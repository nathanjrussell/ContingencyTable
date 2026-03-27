#include "ContingencyTable/FeatureSelector.h"

#include <algorithm>
#include <unordered_set>
#include <stdexcept>

namespace ContingencyTableLib {

FeatureSelector::FeatureSelector() = default;

FeatureSelector::~FeatureSelector() = default;

void FeatureSelector::load(const std::string& path) {
  dataPath_ = path;
  dataTable_.load(dataPath_);
  dataTableLoaded_ = true;

  // Target column is accessed repeatedly during scans; cache it in DataTable.
  dataTable_.setFrequentColumn(static_cast<std::uint64_t>(targetColumn_));
  columnFound_ = false;
  partitionFound_ = false;
}

void FeatureSelector::setTargetColumn(std::size_t columnIndex) {
  targetColumn_ = columnIndex;

  // Hint DataTable that we will repeatedly access the target column.
  if (dataTableLoaded_) {
    if (targetColumn_ >= dataTable_.getColumnCount()) {
      throw std::out_of_range("Target column index is out of range.");
    }
    dataTable_.setFrequentColumn(static_cast<std::uint64_t>(targetColumn_));
  }

  columnFound_ = false;
  partitionFound_ = false;
}

void FeatureSelector::setColumnAlpha(double alpha, bool applyBonferroni) {
  columnAlpha_ = alpha;
  columnBonferroni_ = applyBonferroni;
  columnFound_ = false;
  partitionFound_ = false;
}

void FeatureSelector::setPartitionAlpha(double alpha, bool applyBonferroni) {
  partitionAlpha_ = alpha;
  partitionBonferroni_ = applyBonferroni;
  partitionFound_ = false;
}

void FeatureSelector::setSkipEmptyValues(bool skip) {
  skipEmptyValues_ = skip;
  columnFound_ = false;
  partitionFound_ = false;
}

void FeatureSelector::enabledRows(const std::uint32_t* bitmask, std::size_t sizeInBits) {
  rowBitmask_ = bitmask;
  rowBitmaskSize_ = sizeInBits;
  ownedRowBitmask_.clear();
  columnFound_ = false;
  partitionFound_ = false;
}

void FeatureSelector::enabledRows(
    std::vector<std::uint32_t>& rowIndices,
    std::uint32_t threshold,
    const std::vector<std::uint32_t>& allowedValues,
    std::size_t columnIndex) {
  if (dataPath_.empty() || !dataTableLoaded_) {
    throw std::runtime_error("No data loaded. Call load() first.");
  }

  if (columnIndex >= dataTable_.getColumnCount()) {
    throw std::out_of_range("Column index is out of range.");
  }

  const auto rowCount = static_cast<std::size_t>(dataTable_.getRowCount());
  if (rowIndices.size() < rowCount) {
    throw std::invalid_argument("rowIndices must be at least getRowCount() in size.");
  }

  // Step 1: any marker >= threshold becomes 0.
  for (std::size_t i = 0; i < rowCount; ++i) {
    if (rowIndices[i] >= threshold) {
      rowIndices[i] = 0;
    }
  }

  // Step 2: scan the requested column row-by-row.
  // DataTable convention used by this repo: row 0 is the header row; data starts at row 1.
  // We update rowIndices using the *DataTable row index*.
  // Use a hash set for allowedValues membership to make this O(rows) average instead of O(rows * allowedValues).
  const std::unordered_set<std::uint32_t> allowed(allowedValues.begin(), allowedValues.end());
  for (std::size_t row = 1; row < rowCount; ++row) {
    const auto value = dataTable_.lookupMap(row, static_cast<std::uint64_t>(columnIndex));
    if (allowed.find(value) == allowed.end()) {
      rowIndices[row] = threshold;
    }
  }

  // Step 3: build bitmask (marker==0 => enabled => bit=1).
  const std::size_t sizeInBits = rowCount;
  const std::size_t numWords = (sizeInBits + 31) / 32;

  ownedRowBitmask_.assign(numWords, 0);
  for (std::size_t idx = 0; idx < sizeInBits; ++idx) {
    if (rowIndices[idx] == 0) {
      const std::size_t word = idx / 32;
      const std::size_t bit = idx % 32;
      ownedRowBitmask_[word] |= (1u << bit);
    }
  }

  // Step 4: set this instance's enabled-rows bitmask to the owned storage.
  enabledRows(ownedRowBitmask_.data(), sizeInBits);
}

void FeatureSelector::enabledColumns(const std::uint32_t* bitmask, std::size_t sizeInBits) {
  colBitmask_ = bitmask;
  colBitmaskSize_ = sizeInBits;
  columnFound_ = false;
  partitionFound_ = false;
}

void FeatureSelector::findSignificantColumn() {
  if (dataPath_.empty()) {
    throw std::runtime_error("No data loaded. Call load() first.");
  }
  if (!dataTableLoaded_) {
    throw std::runtime_error("No data loaded. Call load() first.");
  }
  if (colBitmask_ == nullptr) {
    throw std::runtime_error("Column mask not set. Call enabledColumns() first.");
  }

  // Target column is accessed repeatedly during scans; ensure it is pinned in DataTable.
  // (Cheap no-op if already set.)
  dataTable_.setFrequentColumn(static_cast<std::uint64_t>(targetColumn_));

  // Count enabled columns for Bonferroni correction
  const std::size_t numTests = countSetBits(colBitmask_, colBitmaskSize_);
  if (numTests == 0) {
    columnFound_ = false;
    partitionFound_ = false;
    return;
  }

  const double effectiveAlpha = computeEffectiveAlpha(columnAlpha_, columnBonferroni_, numTests);

  // Reset results
  columnFound_ = false;
  partitionFound_ = false;
  bestPValue_ = 1.0;
  bestChiSquare_ = 0.0;
  bestDegreesOfFreedom_ = 0;
  bestColumn_ = 0;

  // Reuse a single loaded ContingencyTable (DataTable::load() is relatively expensive for large datasets).
  ContingencyTable ct;
  ct.load(dataPath_);
  ct.setFirstColumn(targetColumn_);
  ct.setSkipEmptyValues(skipEmptyValues_);
  if (rowBitmask_ != nullptr) {
    ct.setRowFilter(rowBitmask_, rowBitmaskSize_);
  }

  // Test each enabled column
  for (std::size_t col = 0; col < colBitmaskSize_; ++col) {
    if (!isBitSet(colBitmask_, col)) {
      continue;
    }

    ct.setSecondColumn(col);
    ct.build();

    const double pval = ct.getPValue();

    // Track best column (lowest p-value)
    if (pval < bestPValue_) {
      bestPValue_ = pval;
      bestChiSquare_ = ct.getTestStatistic();
      bestDegreesOfFreedom_ = ct.getDegreesOfFreedom();
      bestColumn_ = col;
    }
  }

  // Check if best column is significant
  if (bestPValue_ < effectiveAlpha) {
    columnFound_ = true;

    // Now search for partition on the best column (reuse the already-loaded table instance)
    ct.setSecondColumn(bestColumn_);
    ct.build();

    // Find partition with its own alpha (typically no Bonferroni here)
    const double partitionEffectiveAlpha = computeEffectiveAlpha(partitionAlpha_, partitionBonferroni_, 1);
    auto partition = ct.findOptimalPartition(partitionEffectiveAlpha);

    if (partition.has_value()) {
      partitionFound_ = true;
      partition0_ = std::move(partition->partition0);
      partition1_ = std::move(partition->partition1);
      partitionChiSquare_ = partition->chiSquare;
      partitionPValue_ = partition->pValue;
      partitionDegreesOfFreedom_ = partition->degreesOfFreedom;
    }
  }
}

std::uint64_t FeatureSelector::getRowCount() const {
  if (dataPath_.empty() || !dataTableLoaded_) {
    throw std::runtime_error("No data loaded. Call load() first.");
  }
  return dataTable_.getRowCount();
}

std::uint64_t FeatureSelector::getColumnCount() const {
  if (dataPath_.empty() || !dataTableLoaded_) {
    throw std::runtime_error("No data loaded. Call load() first.");
  }
  return dataTable_.getColumnCount();
}

std::string FeatureSelector::getColumnHeader(std::uint64_t col) const {
  if (dataPath_.empty() || !dataTableLoaded_) {
    throw std::runtime_error("No data loaded. Call load() first.");
  }
  return dataTable_.getColumnHeader(col);
}

bool FeatureSelector::significantColumnFound() const {
  return columnFound_;
}

std::size_t FeatureSelector::getSignificantColumnIndex() const {
  ensureColumnFound();
  return bestColumn_;
}

double FeatureSelector::getColumnPValue() const {
  ensureColumnFound();
  return bestPValue_;
}

std::uint64_t FeatureSelector::getColumnDegreesOfFreedom() const {
  ensureColumnFound();
  return bestDegreesOfFreedom_;
}

double FeatureSelector::getColumnTestStatistic() const {
  ensureColumnFound();
  return bestChiSquare_;
}

bool FeatureSelector::significantPartitionFound() const {
  return partitionFound_;
}

double FeatureSelector::getPartitionPValue() const {
  ensurePartitionFound();
  return partitionPValue_;
}

std::uint64_t FeatureSelector::getPartitionDegreesOfFreedom() const {
  ensurePartitionFound();
  return partitionDegreesOfFreedom_;
}

double FeatureSelector::getPartitionTestStatistic() const {
  ensurePartitionFound();
  return partitionChiSquare_;
}

std::vector<std::uint32_t> FeatureSelector::getFirstPartition() const {
  ensurePartitionFound();
  return partition0_;
}

std::vector<std::uint32_t> FeatureSelector::getSecondPartition() const {
  ensurePartitionFound();
  return partition1_;
}

double FeatureSelector::computeEffectiveAlpha(double alpha, bool applyBonferroni, std::size_t numTests) const {
  if (applyBonferroni && numTests > 1) {
    return alpha / static_cast<double>(numTests);
  }
  return alpha;
}

void FeatureSelector::ensureColumnFound() const {
  if (!columnFound_) {
    throw std::runtime_error("No significant column found. Check significantColumnFound() before calling this method.");
  }
}

void FeatureSelector::ensurePartitionFound() const {
  if (!partitionFound_) {
    throw std::runtime_error("No significant partition found. Check significantPartitionFound() before calling this method.");
  }
}

}  // namespace ContingencyTableLib

