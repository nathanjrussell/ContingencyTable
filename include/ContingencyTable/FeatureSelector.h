#pragma once

#include <ContingencyTable/ContingencyTable.h>
#include <DataTable/DataTable.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace ContingencyTableLib {

class FeatureSelector {
 public:
  FeatureSelector();
  ~FeatureSelector();

  // Configuration
  void load(const std::string& path);
  void setTargetColumn(std::size_t columnIndex);
  void setColumnAlpha(double alpha, bool applyBonferroni);
  void setPartitionAlpha(double alpha, bool applyBonferroni);
  void setSkipEmptyValues(bool skip);

  // Row/column filtering via bitmask
  void enabledRows(const std::uint32_t* bitmask, std::size_t sizeInBits);

  // Overload: derive enabled row bitmask from row-index markers and a column whitelist.
  // - rowIndices: per-row markers (will be mutated). 0 means enabled.
  // - threshold: any marker >= threshold is reset to 0 before applying the column whitelist.
  // - allowedValues: values considered "in set" for the column check.
  // - columnIndex: DataTable column to inspect (row-by-row).
  void enabledRows(
      std::vector<std::uint32_t>& rowIndices,
      std::uint32_t threshold,
      const std::vector<std::uint32_t>& allowedValues,
      std::size_t columnIndex);

  void enabledColumns(const std::uint32_t* bitmask, std::size_t sizeInBits);

  // Execute feature selection
  void findSignificantColumn();

  // DataTable access methods
  std::uint64_t getRowCount() const;
  std::uint64_t getColumnCount() const;
  std::string getColumnHeader(std::uint64_t col) const;

  // Column selection results (throws if !significantColumnFound())
  bool significantColumnFound() const;
  std::size_t getSignificantColumnIndex() const;
  double getColumnPValue() const;
  std::uint64_t getColumnDegreesOfFreedom() const;
  double getColumnTestStatistic() const;

  // Partition search (throws if !significantPartitionFound())
  bool significantPartitionFound() const;
  double getPartitionPValue() const;
  std::uint64_t getPartitionDegreesOfFreedom() const;
  double getPartitionTestStatistic() const;
  std::vector<std::uint32_t> getFirstPartition() const;
  std::vector<std::uint32_t> getSecondPartition() const;

  // Number of included (active) rows whose best-feature value falls into each partition.
  // Partition "one" corresponds to getFirstPartition() and partition "two" to getSecondPartition().
  std::uint64_t getPartitionOneRowCount() const;
  std::uint64_t getPartitionTwoRowCount() const;

 private:
  std::string dataPath_;
    // Cached DataTable (loaded once in load()). Marked mutable so const getters can reuse it.
    mutable DataTableLib::DataTable dataTable_;
    bool dataTableLoaded_ = false;
  std::size_t targetColumn_ = 0;

  // Alpha thresholds
  double columnAlpha_ = 0.05;
  bool columnBonferroni_ = true;
  double partitionAlpha_ = 0.05;
  bool partitionBonferroni_ = false;

  // Row/column filters (raw bitmasks - not owned)
  const std::uint32_t* rowBitmask_ = nullptr;
  std::size_t rowBitmaskSize_ = 0;
  const std::uint32_t* colBitmask_ = nullptr;
  std::size_t colBitmaskSize_ = 0;

  // Owned row bitmask storage for the overload that constructs a mask internally.
  // When non-empty, rowBitmask_ points to this buffer.
  std::vector<std::uint32_t> ownedRowBitmask_;

  bool skipEmptyValues_ = true;

  // Column selection results
  bool columnFound_ = false;
  std::size_t bestColumn_ = 0;
  double bestPValue_ = 1.0;
  double bestChiSquare_ = 0.0;
  std::uint64_t bestDegreesOfFreedom_ = 0;

  // Partition results
  bool partitionFound_ = false;
  double partitionPValue_ = 1.0;
  double partitionChiSquare_ = 0.0;
  std::uint64_t partitionDegreesOfFreedom_ = 0;
  std::vector<std::uint32_t> partition0_;
  std::vector<std::uint32_t> partition1_;

  // Active-row counts for each partition side (computed when partitionFound_ is true).
  std::uint64_t partition0RowCount_ = 0;
  std::uint64_t partition1RowCount_ = 0;

  // Helpers
  double computeEffectiveAlpha(double alpha, bool applyBonferroni, std::size_t numTests) const;
  void ensureColumnFound() const;
  void ensurePartitionFound() const;
};

}  // namespace ContingencyTableLib

