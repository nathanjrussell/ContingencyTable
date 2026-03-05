#include "ContingencyTable/FeatureSelector.h"

#include <stdexcept>

namespace ContingencyTableLib {

FeatureSelector::FeatureSelector() = default;

FeatureSelector::~FeatureSelector() = default;

void FeatureSelector::load(const std::string& path) {
  dataPath_ = path;
  columnFound_ = false;
  partitionFound_ = false;
}

void FeatureSelector::setTargetColumn(std::size_t columnIndex) {
  targetColumn_ = columnIndex;
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
  columnFound_ = false;
  partitionFound_ = false;
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
  if (colBitmask_ == nullptr) {
    throw std::runtime_error("Column mask not set. Call enabledColumns() first.");
  }

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

  // Test each enabled column
  for (std::size_t col = 0; col < colBitmaskSize_; ++col) {
    if (!isBitSet(colBitmask_, col)) {
      continue;
    }

    // Build contingency table for this column vs target
    ContingencyTable ct;
    ct.load(dataPath_);
    ct.setFirstColumn(targetColumn_);
    ct.setSecondColumn(col);
    ct.setSkipEmptyValues(skipEmptyValues_);

    if (rowBitmask_ != nullptr) {
      ct.setRowFilter(rowBitmask_, rowBitmaskSize_);
    }

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

    // Now search for partition on the best column
    ContingencyTable bestTable;
    bestTable.load(dataPath_);
    bestTable.setFirstColumn(targetColumn_);
    bestTable.setSecondColumn(bestColumn_);
    bestTable.setSkipEmptyValues(skipEmptyValues_);

    if (rowBitmask_ != nullptr) {
      bestTable.setRowFilter(rowBitmask_, rowBitmaskSize_);
    }

    bestTable.build();

    // Find partition with its own alpha (typically no Bonferroni here)
    const double partitionEffectiveAlpha = computeEffectiveAlpha(partitionAlpha_, partitionBonferroni_, 1);
    auto partition = bestTable.findOptimalPartition(partitionEffectiveAlpha);

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
  if (dataPath_.empty()) {
    throw std::runtime_error("No data loaded. Call load() first.");
  }
  DataTableLib::DataTable dt;
  dt.load(dataPath_);
  return dt.getRowCount();
}

std::uint64_t FeatureSelector::getColumnCount() const {
  if (dataPath_.empty()) {
    throw std::runtime_error("No data loaded. Call load() first.");
  }
  DataTableLib::DataTable dt;
  dt.load(dataPath_);
  return dt.getColumnCount();
}

std::string FeatureSelector::getColumnHeader(std::uint64_t col) const {
  if (dataPath_.empty()) {
    throw std::runtime_error("No data loaded. Call load() first.");
  }
  DataTableLib::DataTable dt;
  dt.load(dataPath_);
  return dt.getColumnHeader(col);
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

