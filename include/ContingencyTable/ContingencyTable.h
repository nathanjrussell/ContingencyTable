#pragma once

#include <DataTable/DataTable.h>

#include <cstdint>
#include <map>
#include <optional>
#include <utility>
#include <string>
#include <vector>

namespace ContingencyTableLib {

// Reusable bitmask helpers
inline bool isBitSet(const std::uint32_t* bitmask, std::size_t index) {
  const std::size_t word = index / 32;
  const std::size_t bit = index % 32;
  return (bitmask[word] & (1u << bit)) != 0;
}

inline std::size_t countSetBits(const std::uint32_t* bitmask, std::size_t sizeInBits) {
  std::size_t count = 0;
  const std::size_t numWords = (sizeInBits + 31) / 32;
  for (std::size_t i = 0; i < numWords; ++i) {
    count += __builtin_popcount(bitmask[i]);
  }
  return count;
}

struct PartitionResult {
  std::vector<std::uint32_t> partition0;
  std::vector<std::uint32_t> partition1;
  double chiSquare;
  double pValue;
  std::uint64_t degreesOfFreedom;
};

// ...existing code...

class ContingencyTable : public DataTableLib::DataTable {
 public:
  ContingencyTable();
  ~ContingencyTable();

  // Column selection
  void setFirstColumn(std::uint64_t columnId);
  void setFirstColumn(const std::string& header);
  void setSecondColumn(std::uint64_t columnId);
  void setSecondColumn(const std::string& header);

  // Row filtering via bitmask
  void setRowFilter(const std::uint32_t* bitmask, std::size_t sizeInBits);

  // Empty value policy (default: true = skip rows with 0 mapped values)
  void setSkipEmptyValues(bool skip);

  // Build the contingency table
  void build();

  // Results
  double getTestStatistic() const;
  std::uint64_t getDegreesOfFreedom() const;
  double getPValue() const;

  // Find optimal partition of column B (requires build() to be called first)
  std::optional<PartitionResult> findOptimalPartition(double alpha);

  // Partition results (available after successful findOptimalPartition() call)
  double getPartitionChiSquare() const;
  double getPartitionPValue() const;
  std::uint64_t getPartitionDegreesOfFreedom() const;

  // Count how many included (active) rows fall into each side of a partition of column B.
  // Requires build() to have been called so jointCounts_ is populated.
  // Returns {rowsInPartition0, rowsInPartition1}.
  std::pair<std::uint64_t, std::uint64_t> countRowsInPartitions(
      const std::vector<std::uint32_t>& partition0) const;

 private:
  std::uint64_t colA_ = 0;
  std::uint64_t colB_ = 0;
  const std::uint32_t* rowBitmask_ = nullptr;
  std::size_t bitmaskSize_ = 0;
  bool skipEmptyValues_ = true;
  bool dirty_ = true;

  double chiSquare_ = 0.0;
  std::uint64_t degreesOfFreedom_ = 0;

  // Partition results
  double partitionChiSquare_ = 0.0;
  double partitionPValue_ = 1.0;
  std::uint64_t partitionDegreesOfFreedom_ = 0;

  // Joint distribution table: (featureA, featureB) -> count
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint64_t> jointCounts_;

  void invalidate();
  bool isRowIncluded(std::uint64_t row) const;
  double chiSquarePValue(double chiSq, std::uint64_t df) const;
  void buildJointCounts();

  // Partition search helpers
  std::uint32_t determineRestarts(std::size_t numFeaturesB, double adjustedAlpha) const;
  std::optional<PartitionResult> greedySearch(double adjustedAlpha);
};

}  // namespace ContingencyTableLib

