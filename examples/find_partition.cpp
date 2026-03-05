#include "ContingencyTable/ContingencyTable.h"

#include <iostream>
#include <vector>
#include <iomanip>

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "Usage: find_partition <datatable_dir> <col1> <col2> [alpha=0.05]\n";
    return 1;
  }

  try {
    const std::string outputDir = argv[1];
    const auto col1 = static_cast<std::uint64_t>(std::strtoull(argv[2], nullptr, 10));
    const auto col2 = static_cast<std::uint64_t>(std::strtoull(argv[3], nullptr, 10));
    const double alpha = (argc >= 5) ? std::atof(argv[4]) : 0.05;

    ContingencyTableLib::ContingencyTable table;
    table.load(outputDir);

    std::cout << "Loaded DataTable: " << table.getRowCount() << " rows, "
              << table.getColumnCount() << " columns\n";

    table.setFirstColumn(col1);
    table.setSecondColumn(col2);
    table.setSkipEmptyValues(true);
    table.build();

    std::cout << "\n=== ORIGINAL CONTINGENCY TABLE ===\n";
    std::cout << "Chi-square: " << std::fixed << std::setprecision(6)
              << table.getTestStatistic() << "\n";
    std::cout << "P-value:    " << table.getPValue() << "\n";
    std::cout << "DF:         " << table.getDegreesOfFreedom() << "\n";

    std::cout << "\n=== SEARCHING FOR OPTIMAL PARTITION ===\n";
    auto partition = table.findOptimalPartition(alpha);

    if (partition) {
      std::cout << "Partition FOUND and is significant!\n";
      std::cout << "\nPartition 0 features: ";
      for (auto f : partition->partition0) {
        std::cout << f << " ";
      }
      std::cout << "\nPartition 1 features: ";
      for (auto f : partition->partition1) {
        std::cout << f << " ";
      }
      std::cout << "\n=== PARTITION CONTINGENCY TABLE ===\n";
      std::cout << "Chi-square: " << std::fixed << std::setprecision(6)
                << table.getPartitionChiSquare() << "\n";
      std::cout << "P-value:    " << table.getPartitionPValue() << "\n";
      std::cout << "DF:         " << table.getPartitionDegreesOfFreedom() << "\n";

      // Debug: Print which features are in each partition
      std::cout << "\nDebug - Partition contents:\n";
      std::cout << "  Partition 0 (feature values): ";
      for (auto f : partition->partition0) {
        std::cout << "'" << table.getColumnValue(col2, f) << "' ";
      }
      std::cout << "\n  Partition 1 (feature values): ";
      for (auto f : partition->partition1) {
        std::cout << "'" << table.getColumnValue(col2, f) << "' ";
      }
      std::cout << "\n";

      // Manually compute chi-square to verify
      std::cout << "\nManual verification of chi-square:\n";
      std::cout << "  Partition 0 features: ";
      for (auto f : partition->partition0) std::cout << f << " ";
      std::cout << "\n  Partition 1 features: ";
      for (auto f : partition->partition1) std::cout << f << " ";
      std::cout << "\n";

      // Output JSON format for Python parsing
      std::cout << "\n=== JSON OUTPUT ===\n";
      std::cout << "{\n";
      std::cout << "  \"partition0\": [";
      for (size_t i = 0; i < partition->partition0.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << partition->partition0[i];
      }
      std::cout << "],\n";
      std::cout << "  \"partition1\": [";
      for (size_t i = 0; i < partition->partition1.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << partition->partition1[i];
      }
      std::cout << "],\n";
      std::cout << "  \"chi_square\": " << std::fixed << std::setprecision(6)
                << table.getPartitionChiSquare() << ",\n";
      std::cout << "  \"p_value\": " << std::fixed << std::setprecision(9)
                << table.getPartitionPValue() << ",\n";
      std::cout << "  \"df\": " << table.getPartitionDegreesOfFreedom() << "\n";
      std::cout << "}\n";
    } else {
      std::cout << "No significant partition found at alpha = " << alpha << "\n";
    }

    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 2;
  }
}



