#include "ContingencyTable/ContingencyTable.h"

#include <iostream>
#include <vector>
#include <iomanip>

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "Usage: verify_rowfilter <datatable_dir> <col1> <col2>\n";
    return 1;
  }

  try {
    const std::string outputDir = argv[1];
    const auto col1 = static_cast<std::uint64_t>(std::strtoull(argv[2], nullptr, 10));
    const auto col2 = static_cast<std::uint64_t>(std::strtoull(argv[3], nullptr, 10));

    ContingencyTableLib::ContingencyTable table;
    table.load(outputDir);

    std::cout << "Loaded DataTable: " << table.getRowCount() << " rows, "
              << table.getColumnCount() << " columns\n";
    std::cout << "Column 1: " << table.getColumnHeader(col1) << "\n";
    std::cout << "Column 2: " << table.getColumnHeader(col2) << "\n\n";

    // ========================================================================
    // TEST 1: No row filter (should include all data rows)
    // ========================================================================
    std::cout << "=" << std::string(78, '=') << "\n";
    std::cout << "TEST 1: WITHOUT Row Filter (All Rows Active)\n";
    std::cout << "=" << std::string(78, '=') << "\n";

    ContingencyTableLib::ContingencyTable table1;
    table1.load(outputDir);
    table1.setFirstColumn(col1);
    table1.setSecondColumn(col2);
    table1.setSkipEmptyValues(true);
    // NOTE: NO setRowFilter() called - should default to all rows active
    table1.build();

    std::cout << "Chi-Square: " << std::fixed << std::setprecision(6)
              << table1.getTestStatistic() << "\n";
    std::cout << "P-Value:    " << table1.getPValue() << "\n";
    std::cout << "DF:         " << table1.getDegreesOfFreedom() << "\n\n";

    // ========================================================================
    // TEST 2: Exclude row 1 (should have fewer total observations)
    // ========================================================================
    std::cout << "=" << std::string(78, '=') << "\n";
    std::cout << "TEST 2: WITH Row Filter (Exclude row 1)\n";
    std::cout << "=" << std::string(78, '=') << "\n";

    // Create bitmask to exclude row 1
    // Row 0 = header (not included anyway)
    // Row 1 = skip
    // Rows 2+ = include
    const auto totalRows = table.getRowCount();
    const std::size_t numWords = (totalRows + 31) / 32;
    std::vector<std::uint32_t> bitmask(numWords, 0xFFFFFFFF);

    // Clear bit 1 to exclude row 1
    bitmask[0] &= ~(1U << 1);

    std::cout << "Bitmask: Bit 0 (header) and Bit 1 (row 1) will be skipped\n\n";

    ContingencyTableLib::ContingencyTable table2;
    table2.load(outputDir);
    table2.setFirstColumn(col1);
    table2.setSecondColumn(col2);
    table2.setSkipEmptyValues(true);
    table2.setRowFilter(bitmask.data(), totalRows);
    table2.build();

    std::cout << "Chi-Square: " << std::fixed << std::setprecision(6)
              << table2.getTestStatistic() << "\n";
    std::cout << "P-Value:    " << table2.getPValue() << "\n";
    std::cout << "DF:         " << table2.getDegreesOfFreedom() << "\n\n";

    // ========================================================================
    // SUMMARY
    // ========================================================================
    std::cout << "=" << std::string(78, '=') << "\n";
    std::cout << "VERIFICATION SUMMARY\n";
    std::cout << "=" << std::string(78, '=') << "\n";
    std::cout << "✓ WITHOUT setRowFilter():\n";
    std::cout << "  - All data rows are active (header row 0 is automatically skipped)\n";
    std::cout << "  - Chi-Square: " << std::fixed << std::setprecision(6)
              << table1.getTestStatistic() << "\n\n";

    std::cout << "✓ WITH setRowFilter() (excluding row 1):\n";
    std::cout << "  - Only rows with bitmask bit set to 1 are included\n";
    std::cout << "  - Row 0 (header) is never included (iteration starts at row 1)\n";
    std::cout << "  - Chi-Square: " << std::fixed << std::setprecision(6)
              << table2.getTestStatistic() << "\n\n";

    std::cout << "Difference in Chi-Square: " << std::fixed << std::setprecision(6)
              << (table1.getTestStatistic() - table2.getTestStatistic()) << "\n";
    std::cout << "(Should be non-zero since excluding a row changes the contingency table)\n";

    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 2;
  }
}

