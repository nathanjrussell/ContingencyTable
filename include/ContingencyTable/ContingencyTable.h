#pragma once

#include <DataTable/DataTable.h>

#include <cstdint>
#include <string>

namespace ContingencyTableLib {

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

 private:
  std::uint64_t colA_ = 0;
  std::uint64_t colB_ = 0;
  const std::uint32_t* rowBitmask_ = nullptr;
  std::size_t bitmaskSize_ = 0;
  bool skipEmptyValues_ = true;
  bool dirty_ = true;

  double chiSquare_ = 0.0;
  std::uint64_t degreesOfFreedom_ = 0;

  void invalidate();
  bool isRowIncluded(std::uint64_t row) const;
  double chiSquarePValue(double chiSq, std::uint64_t df) const;
};

}  // namespace ContingencyTableLib

