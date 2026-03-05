#include "ContingencyTable/ContingencyTable.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

namespace ContingencyTableLib {

ContingencyTable::ContingencyTable() = default;

ContingencyTable::~ContingencyTable() = default;

void ContingencyTable::setFirstColumn(std::uint64_t columnId) {
  if (columnId >= getColumnCount()) {
    throw std::out_of_range("First column ID is out of range.");
  }
  colA_ = columnId;
  invalidate();
}

void ContingencyTable::setFirstColumn(const std::string& header) {
  colA_ = getColumnIndex(header);
  invalidate();
}

void ContingencyTable::setSecondColumn(std::uint64_t columnId) {
  if (columnId >= getColumnCount()) {
    throw std::out_of_range("Second column ID is out of range.");
  }
  colB_ = columnId;
  invalidate();
}

void ContingencyTable::setSecondColumn(const std::string& header) {
  colB_ = getColumnIndex(header);
  invalidate();
}

void ContingencyTable::setRowFilter(const std::uint32_t* bitmask, std::size_t sizeInBits) {
  rowBitmask_ = bitmask;
  bitmaskSize_ = sizeInBits;
  invalidate();
}

void ContingencyTable::setSkipEmptyValues(bool skip) {
  skipEmptyValues_ = skip;
  invalidate();
}

void ContingencyTable::invalidate() {
  dirty_ = true;
  chiSquare_ = 0.0;
  degreesOfFreedom_ = 0;
}

bool ContingencyTable::isRowIncluded(std::uint64_t row) const {
  if (rowBitmask_ == nullptr) {
    return true;
  }
  if (row >= bitmaskSize_) {
    return false;
  }
  const std::uint32_t wordIndex = static_cast<std::uint32_t>(row / 32);
  const std::uint32_t bitIndex = static_cast<std::uint32_t>(row % 32);
  return (rowBitmask_[wordIndex] & (1U << bitIndex)) != 0;
}

void ContingencyTable::build() {
  chiSquare_ = 0.0;
  degreesOfFreedom_ = 0;

  const auto rowCount = getRowCount();
  if (rowCount <= 1) {
    dirty_ = false;
    return;
  }

  // Count observed (featureA, featureB) pairs
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint64_t> counts;
  std::map<std::uint32_t, std::uint32_t> featureMapA;
  std::map<std::uint32_t, std::uint32_t> featureMapB;

  for (std::uint64_t row = 1; row < rowCount; ++row) {
    if (!isRowIncluded(row)) {
      continue;
    }

    const auto featureA = lookupMap(row, colA_);
    const auto featureB = lookupMap(row, colB_);

    if (skipEmptyValues_ && (featureA == 0 || featureB == 0)) {
      continue;
    }

    counts[{featureA, featureB}] += 1;
  }

  if (counts.empty()) {
    dirty_ = false;
    return;
  }

  // Build dense remapping for observed features only
  std::uint32_t nextA = 0;
  std::uint32_t nextB = 0;

  for (const auto& kv : counts) {
    const auto featureA = kv.first.first;
    const auto featureB = kv.first.second;

    if (featureMapA.find(featureA) == featureMapA.end()) {
      featureMapA[featureA] = nextA++;
    }
    if (featureMapB.find(featureB) == featureMapB.end()) {
      featureMapB[featureB] = nextB++;
    }
  }

  const std::size_t numRowCategories = featureMapA.size();
  const std::size_t numColCategories = featureMapB.size();

  if (numRowCategories < 2 || numColCategories < 2) {
    dirty_ = false;
    return;
  }

  // Build dense observed table
  std::vector<std::vector<std::uint64_t>> observed(
      numRowCategories,
      std::vector<std::uint64_t>(numColCategories, 0));

  for (const auto& kv : counts) {
    const auto denseA = featureMapA.at(kv.first.first);
    const auto denseB = featureMapB.at(kv.first.second);
    observed[denseA][denseB] = kv.second;
  }

  // Compute row/col totals
  std::vector<double> rowTotals(numRowCategories, 0.0);
  std::vector<double> colTotals(numColCategories, 0.0);
  double grandTotal = 0.0;

  for (std::size_t r = 0; r < numRowCategories; ++r) {
    for (std::size_t c = 0; c < numColCategories; ++c) {
      const double val = static_cast<double>(observed[r][c]);
      rowTotals[r] += val;
      colTotals[c] += val;
      grandTotal += val;
    }
  }

  if (grandTotal == 0.0) {
    dirty_ = false;
    return;
  }

  // Compute Pearson chi-square statistic
  chiSquare_ = 0.0;
  for (std::size_t r = 0; r < numRowCategories; ++r) {
    for (std::size_t c = 0; c < numColCategories; ++c) {
      const double expected = (rowTotals[r] * colTotals[c]) / grandTotal;
      if (expected > 0.0) {
        const double obs = static_cast<double>(observed[r][c]);
        const double diff = obs - expected;
        chiSquare_ += (diff * diff) / expected;
      }
    }
  }

  degreesOfFreedom_ = (numRowCategories - 1) * (numColCategories - 1);
  dirty_ = false;
}

double ContingencyTable::getTestStatistic() const {
  if (dirty_) {
    throw std::logic_error("Must call build() before reading test statistic.");
  }
  return chiSquare_;
}

std::uint64_t ContingencyTable::getDegreesOfFreedom() const {
  if (dirty_) {
    throw std::logic_error("Must call build() before reading degrees of freedom.");
  }
  return degreesOfFreedom_;
}

double ContingencyTable::getPValue() const {
  if (dirty_) {
    throw std::logic_error("Must call build() before reading p-value.");
  }
  return chiSquarePValue(chiSquare_, degreesOfFreedom_);
}

// Built-in chi-square p-value approximation (Wilson-Hilferty + normal tail)
double ContingencyTable::chiSquarePValue(double chiSq, std::uint64_t df) const {
  if (df == 0 || chiSq < 0.0) {
    return 1.0;
  }

  const double dfDouble = static_cast<double>(df);

  // Wilson-Hilferty transform to approximate normality
  const double z = std::pow(chiSq / dfDouble, 1.0 / 3.0)
                 - (1.0 - 2.0 / (9.0 * dfDouble));
  const double zNorm = z / std::sqrt(2.0 / (9.0 * dfDouble));

  // Upper tail probability via complementary error function
  // P(Z > zNorm) = 0.5 * erfc(zNorm / sqrt(2))
  const double pValue = 0.5 * std::erfc(zNorm / std::sqrt(2.0));

  return std::max(0.0, std::min(1.0, pValue));
}

}  // namespace ContingencyTableLib
