#include "ContingencyTable/ContingencyTable.h"

#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "Usage: contingency_example <datatable_dir> <col1_index> <col2_index>\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  ./contingency_example titanic_output 1 3\n";
    return 1;
  }

  try {
    const std::string outputDir = argv[1];
    const auto col1 = static_cast<std::uint64_t>(std::strtoull(argv[2], nullptr, 10));
    const auto col2 = static_cast<std::uint64_t>(std::strtoull(argv[3], nullptr, 10));

    // Load the parsed DataTable
    ContingencyTableLib::ContingencyTable table;
    table.load(outputDir);

    std::cout << "Loaded DataTable: " << table.getRowCount() << " rows, "
              << table.getColumnCount() << " columns\n\n";

    // Set up the analysis
    table.setFirstColumn(col1);
    table.setSecondColumn(col2);
    table.setSkipEmptyValues(true);  // Skip rows with missing values

    // Build the contingency table
    table.build();

    // Display results
    std::cout << "Chi-square statistic: " << table.getTestStatistic() << "\n";
    std::cout << "Degrees of freedom:   " << table.getDegreesOfFreedom() << "\n";
    std::cout << "P-value:              " << table.getPValue() << "\n";

    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 2;
  }
}

