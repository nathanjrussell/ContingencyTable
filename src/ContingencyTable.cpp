#include "ContingencyTable/ContingencyTable.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <unordered_set>
#include <utility>
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
  return isBitSet(rowBitmask_, row);
}

void ContingencyTable::build() {
  chiSquare_ = 0.0;
  degreesOfFreedom_ = 0;
  jointCounts_.clear();

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

  // Store joint counts for later use in findOptimalPartition
  jointCounts_ = counts;

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

double ContingencyTable::getPartitionChiSquare() const {
  return partitionChiSquare_;
}

double ContingencyTable::getPartitionPValue() const {
  return partitionPValue_;
}

std::uint64_t ContingencyTable::getPartitionDegreesOfFreedom() const {
  return partitionDegreesOfFreedom_;
}

std::pair<std::uint64_t, std::uint64_t> ContingencyTable::countRowsInPartitions(
    const std::vector<std::uint32_t>& partition0) const {
  if (dirty_) {
    throw std::logic_error("Must call build() before counting partition row totals.");
  }

  if (jointCounts_.empty()) {
    return {0ULL, 0ULL};
  }

  const std::unordered_set<std::uint32_t> p0(partition0.begin(), partition0.end());

  std::uint64_t count0 = 0;
  std::uint64_t count1 = 0;

  for (const auto& kv : jointCounts_) {
    const std::uint32_t featureB = kv.first.second;
    const std::uint64_t n = kv.second;

    if (p0.find(featureB) != p0.end()) {
      count0 += n;
    } else {
      count1 += n;
    }
  }

  return {count0, count1};
}

std::map<std::uint32_t, std::uint64_t> ContingencyTable::getFirstColumnCounts() const {
  if (dirty_) {
    throw std::logic_error("Must call build() before reading first-column counts.");
  }

  std::map<std::uint32_t, std::uint64_t> countsA;
  for (const auto& kv : jointCounts_) {
    const std::uint32_t featureA = kv.first.first;
    const std::uint64_t n = kv.second;
    countsA[featureA] += n;
  }

  return countsA;
}

std::uint32_t ContingencyTable::determineRestarts(std::size_t numFeaturesB, double adjustedAlpha) const {
  // Bonferroni correction: assume testing ~2^|B| partitions
  // Map adjusted alpha (post-correction) to restart count

  std::uint32_t restarts;

  if (adjustedAlpha < 0.0001) {
    restarts = 15;
  } else if (adjustedAlpha < 0.001) {
    restarts = 10;
  } else if (adjustedAlpha < 0.01) {
    restarts = 7;
  } else {
    restarts = 4;
  }

  // Cap at 5 for very large feature sets (cost scaling)
  if (numFeaturesB > 30) {
    restarts = std::min(restarts, 5U);
  }

  return restarts;
}

std::optional<PartitionResult> ContingencyTable::greedySearch(double adjustedAlpha) {
  if (jointCounts_.empty()) {
    return std::nullopt;
  }

  // Extract unique features in column B
  std::set<std::uint32_t> featuresB;
  for (const auto& kv : jointCounts_) {
    featuresB.insert(kv.first.second);
  }

  if (featuresB.size() < 2) {
    return std::nullopt;
  }

  // Initialize: start with random partition
  std::vector<std::uint32_t> partition0(featuresB.begin(), featuresB.end());
  std::vector<std::uint32_t> partition1;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(partition0.begin(), partition0.end(), gen);

  // Move half to partition1
  std::size_t midpoint = partition0.size() / 2;
  partition1.assign(partition0.begin() + midpoint, partition0.end());
  partition0.erase(partition0.begin() + midpoint, partition0.end());

  PartitionResult bestResult;
  bestResult.partition0 = partition0;
  bestResult.partition1 = partition1;

  // Evaluate initial partition
  {
    std::set<std::uint32_t> featuresInP0(partition0.begin(), partition0.end());
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint64_t> repartitioned;

    for (const auto& kv : jointCounts_) {
      auto featureB = kv.first.second;
      auto mappedB = (featuresInP0.count(featureB) > 0) ? 0U : 1U;
      repartitioned[{kv.first.first, mappedB}] += kv.second;
    }

    // Compute chi-square on repartitioned table
    std::map<std::uint32_t, std::uint32_t> featureMapA;
    std::map<std::uint32_t, std::uint32_t> featureMapB;
    std::uint32_t nextA = 0, nextB = 0;

    for (const auto& kv : repartitioned) {
      if (featureMapA.find(kv.first.first) == featureMapA.end()) {
        featureMapA[kv.first.first] = nextA++;
      }
      if (featureMapB.find(kv.first.second) == featureMapB.end()) {
        featureMapB[kv.first.second] = nextB++;
      }
    }

    const std::size_t numRowCategories = featureMapA.size();
    const std::size_t numColCategories = featureMapB.size();

    if (numRowCategories >= 2 && numColCategories >= 2) {
      std::vector<std::vector<std::uint64_t>> observed(
          numRowCategories,
          std::vector<std::uint64_t>(numColCategories, 0));

      for (const auto& kv : repartitioned) {
        const auto denseA = featureMapA.at(kv.first.first);
        const auto denseB = featureMapB.at(kv.first.second);
        observed[denseA][denseB] = kv.second;
      }

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

      double chiSq = 0.0;
      if (grandTotal > 0.0) {
      for (std::size_t r = 0; r < numRowCategories; ++r) {
        for (std::size_t c = 0; c < numColCategories; ++c) {
          const double expected = (rowTotals[r] * colTotals[c]) / grandTotal;
          if (expected > 0.0) {
            const double obs = static_cast<double>(observed[r][c]);
            const double diff = obs - expected;
            chiSq += (diff * diff) / expected;
          }
        }
      }

        bestResult.chiSquare = chiSq;
        const std::uint64_t df = (numRowCategories - 1) * (numColCategories - 1);
        bestResult.pValue = chiSquarePValue(chiSq, df);
        bestResult.degreesOfFreedom = df;
      }
    }
  }

  // Greedy refinement: try moving each feature
  bool improved = true;
  while (improved) {
    improved = false;
    double bestChiSq = bestResult.chiSquare;

    for (std::size_t i = 0; i < partition0.size(); ++i) {
      // Try moving partition0[i] to partition1
      std::uint32_t feature = partition0[i];
      std::set<std::uint32_t> featuresInP0(partition0.begin(), partition0.end());
      featuresInP0.erase(feature);

      std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint64_t> repartitioned;
      for (const auto& kv : jointCounts_) {
        auto featureB = kv.first.second;
        auto mappedB = (featuresInP0.count(featureB) > 0) ? 0U : 1U;
        repartitioned[{kv.first.first, mappedB}] += kv.second;
      }

      // Compute chi-square (similar logic as above)
      std::map<std::uint32_t, std::uint32_t> featureMapA;
      std::map<std::uint32_t, std::uint32_t> featureMapB;
      std::uint32_t nextA = 0, nextB = 0;

      for (const auto& kv : repartitioned) {
        if (featureMapA.find(kv.first.first) == featureMapA.end()) {
          featureMapA[kv.first.first] = nextA++;
        }
        if (featureMapB.find(kv.first.second) == featureMapB.end()) {
          featureMapB[kv.first.second] = nextB++;
        }
      }

      const std::size_t numRowCategories = featureMapA.size();
      const std::size_t numColCategories = featureMapB.size();

      if (numRowCategories >= 2 && numColCategories >= 2) {
        std::vector<std::vector<std::uint64_t>> observed(
            numRowCategories,
            std::vector<std::uint64_t>(numColCategories, 0));

        for (const auto& kv : repartitioned) {
          const auto denseA = featureMapA.at(kv.first.first);
          const auto denseB = featureMapB.at(kv.first.second);
          observed[denseA][denseB] = kv.second;
        }

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

        if (grandTotal > 0.0) {
          double chiSq = 0.0;
          for (std::size_t r = 0; r < numRowCategories; ++r) {
            for (std::size_t c = 0; c < numColCategories; ++c) {
              const double expected = (rowTotals[r] * colTotals[c]) / grandTotal;
              if (expected > 0.0) {
                const double obs = static_cast<double>(observed[r][c]);
                const double diff = obs - expected;
                chiSq += (diff * diff) / expected;
              }
            }
          }

          if (chiSq > bestChiSq) {
            bestChiSq = chiSq;
            partition0.erase(partition0.begin() + i);
            partition1.push_back(feature);
            improved = true;
            const std::uint64_t df = (numRowCategories - 1) * (numColCategories - 1);
            bestResult.chiSquare = chiSq;
            bestResult.pValue = chiSquarePValue(chiSq, df);
            bestResult.degreesOfFreedom = df;
            bestResult.partition0 = partition0;
            bestResult.partition1 = partition1;
            break;
          }
        }
      }
    }

    // Also try moving from partition1 to partition0
    if (!improved) {
      for (std::size_t i = 0; i < partition1.size(); ++i) {
        std::uint32_t feature = partition1[i];
        std::set<std::uint32_t> featuresInP0(partition0.begin(), partition0.end());
        featuresInP0.insert(feature);

        std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint64_t> repartitioned;
        for (const auto& kv : jointCounts_) {
          auto featureB = kv.first.second;
          auto mappedB = (featuresInP0.count(featureB) > 0) ? 0U : 1U;
          repartitioned[{kv.first.first, mappedB}] += kv.second;
        }

        std::map<std::uint32_t, std::uint32_t> featureMapA;
        std::map<std::uint32_t, std::uint32_t> featureMapB;
        std::uint32_t nextA = 0, nextB = 0;

        for (const auto& kv : repartitioned) {
          if (featureMapA.find(kv.first.first) == featureMapA.end()) {
            featureMapA[kv.first.first] = nextA++;
          }
          if (featureMapB.find(kv.first.second) == featureMapB.end()) {
            featureMapB[kv.first.second] = nextB++;
          }
        }

        const std::size_t numRowCategories = featureMapA.size();
        const std::size_t numColCategories = featureMapB.size();

        if (numRowCategories >= 2 && numColCategories >= 2) {
          std::vector<std::vector<std::uint64_t>> observed(
              numRowCategories,
              std::vector<std::uint64_t>(numColCategories, 0));

          for (const auto& kv : repartitioned) {
            const auto denseA = featureMapA.at(kv.first.first);
            const auto denseB = featureMapB.at(kv.first.second);
            observed[denseA][denseB] = kv.second;
          }

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

          if (grandTotal > 0.0) {
            double chiSq = 0.0;
            for (std::size_t r = 0; r < numRowCategories; ++r) {
              for (std::size_t c = 0; c < numColCategories; ++c) {
                const double expected = (rowTotals[r] * colTotals[c]) / grandTotal;
                if (expected > 0.0) {
                  const double obs = static_cast<double>(observed[r][c]);
                  const double diff = obs - expected;
                  chiSq += (diff * diff) / expected;
                }
              }
            }

            if (chiSq > bestChiSq) {
              bestChiSq = chiSq;
              partition1.erase(partition1.begin() + i);
              partition0.push_back(feature);
              improved = true;
              bestResult.chiSquare = chiSq;
              const std::uint64_t df = (numRowCategories - 1) * (numColCategories - 1);
              bestResult.pValue = chiSquarePValue(chiSq, df);
              bestResult.degreesOfFreedom = df;
              bestResult.partition0 = partition0;
              bestResult.partition1 = partition1;
              break;
            }
          }
        }
      }
    }
  }

  return bestResult;
}

std::optional<PartitionResult> ContingencyTable::findOptimalPartition(double alpha) {
  if (dirty_) {
    throw std::logic_error("Must call build() before findOptimalPartition().");
  }

  if (jointCounts_.empty()) {
    return std::nullopt;
  }

  // Extract unique features in column B
  std::set<std::uint32_t> featuresB;
  for (const auto& kv : jointCounts_) {
    featuresB.insert(kv.first.second);
  }

  const std::size_t numFeaturesB = featuresB.size();
  if (numFeaturesB < 2) {
    return std::nullopt;
  }

  // Apply Bonferroni correction: alpha_adjusted = alpha / 2^|B|
  const double adjustedAlpha = alpha / std::pow(2.0, static_cast<double>(numFeaturesB));

  // Determine restart count based on adjusted alpha
  const std::uint32_t restarts = determineRestarts(numFeaturesB, adjustedAlpha);

  // Run greedy search with multiple restarts
  PartitionResult bestResult;
  bestResult.chiSquare = -1.0;
  bestResult.pValue = 1.0;

  for (std::uint32_t i = 0; i < restarts; ++i) {
    auto result = greedySearch(adjustedAlpha);
    if (result && result->chiSquare > bestResult.chiSquare) {
      bestResult = *result;
    }
  }

  // Return result only if it meets the adjusted alpha threshold
  if (bestResult.pValue < adjustedAlpha) {
    partitionChiSquare_ = bestResult.chiSquare;
    partitionPValue_ = bestResult.pValue;
    partitionDegreesOfFreedom_ = bestResult.degreesOfFreedom;
    return bestResult;
  }

  return std::nullopt;
}

}  // namespace ContingencyTableLib
