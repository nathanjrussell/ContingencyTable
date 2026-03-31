#include <ContingencyTable/FeatureSelector.h>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: feature_selector <parsed_datatable_dir>\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  First parse a CSV:\n";
    std::cerr << "    ./parse_csv data.csv parsed_output\n";
    std::cerr << "  Then run feature selection:\n";
    std::cerr << "    ./feature_selector parsed_output\n";
    return 1;
  }

  try {
    const std::string outputDir = argv[1];

    ContingencyTableLib::FeatureSelector fs;
    fs.load(outputDir);

    std::cout << "Loaded DataTable from: " << outputDir << "\n";
    std::cout << "Total rows: " << fs.getRowCount() << "\n";
    std::cout << "Total columns: " << fs.getColumnCount() << "\n\n";

    // ========================================================================
    // Example 1: Simple feature selection (all rows, all columns except 0)
    // ========================================================================
    std::cout << "========== Example 1: Feature Selection ==========\n";
    std::cout << "Target: Column 0 (first column)\n";
    std::cout << "Candidates: All other columns\n\n";

    ContingencyTableLib::FeatureSelector fs1;
    fs1.load(outputDir);
    fs1.setTargetColumn(0);

    // Enable all rows
    const auto rowCount = fs1.getRowCount();
    const std::size_t numWords = (rowCount + 31) / 32;
    std::vector<std::uint32_t> allRowsMask(numWords, 0xFFFFFFFF);
    fs1.enabledRows(allRowsMask.data(), rowCount);

    // Enable all columns except 0 (target column)
    const auto colCount = fs1.getColumnCount();
    const std::size_t colWords = (colCount + 31) / 32;
    std::vector<std::uint32_t> candidateColsMask(colWords, 0xFFFFFFFF);
    candidateColsMask[0] &= ~1U;  // Disable column 0 (the target)
    fs1.enabledColumns(candidateColsMask.data(), colCount);

    // Configure alpha with Bonferroni correction for multiple column testing
    fs1.setColumnAlpha(0.05, true);    // Bonferroni correction for all tested columns
    fs1.setPartitionAlpha(0.05, false); // No correction for partition search
    fs1.setSkipEmptyValues(true);       // Skip rows with missing values

    std::cout << "Searching for most significant feature...\n";
    fs1.findSignificantColumn();

    if (fs1.significantColumnFound()) {
      std::cout << "✓ Significant feature FOUND!\n";
      std::cout << "  Column index: " << fs1.getSignificantColumnIndex() << "\n";
      std::cout << "  Column header: " << fs1.getColumnHeader(fs1.getSignificantColumnIndex()) << "\n";
      std::cout << "  Chi-square: " << fs1.getColumnTestStatistic() << "\n";
      std::cout << "  Degrees of freedom: " << fs1.getColumnDegreesOfFreedom() << "\n";
      std::cout << "  P-value: " << fs1.getColumnPValue() << "\n";

      if (fs1.significantPartitionFound()) {
        std::cout << "\n✓ Optimal partition FOUND!\n";
        auto p0 = fs1.getFirstPartition();
        auto p1 = fs1.getSecondPartition();
        std::cout << "  Partition 0 size: " << p0.size() << " features\n";
        std::cout << "  Partition 1 size: " << p1.size() << " features\n";
        std::cout << "  Partition 0 row count (active rows): " << fs1.getPartitionOneRowCount() << "\n";
        std::cout << "  Partition 1 row count (active rows): " << fs1.getPartitionTwoRowCount() << "\n";
        std::cout << "  Partition chi-square: " << fs1.getPartitionTestStatistic() << "\n";
        std::cout << "  Partition p-value: " << fs1.getPartitionPValue() << "\n";
        std::cout << "  Partition DF: " << fs1.getPartitionDegreesOfFreedom() << "\n";
      } else {
        std::cout << "\n✗ No significant partition found\n";
      }
    } else {
      std::cout << "✗ No significant feature found at alpha=0.05\n";
    }

    // ========================================================================
    // Example 2: Limited feature set (select specific columns)
    // ========================================================================
    std::cout << "\n========== Example 2: Select Specific Columns ==========\n";

    if (colCount >= 4) {
      ContingencyTableLib::FeatureSelector fs2;
      fs2.load(outputDir);
      fs2.setTargetColumn(1);  // Target: column 1

      // Enable all rows
      fs2.enabledRows(allRowsMask.data(), rowCount);

      // Enable only columns 2 and 3 as candidates
      std::vector<std::uint32_t> limitedColsMask(colWords, 0);
      limitedColsMask[0] |= (1U << 2);  // Enable column 2
      limitedColsMask[0] |= (1U << 3);  // Enable column 3
      fs2.enabledColumns(limitedColsMask.data(), colCount);

      std::cout << "Target: Column 1 (" << fs2.getColumnHeader(1) << ")\n";
      std::cout << "Candidates: Column 2 (" << fs2.getColumnHeader(2) << "), "
                << "Column 3 (" << fs2.getColumnHeader(3) << ")\n\n";

      fs2.setColumnAlpha(0.05, true);    // Bonferroni for 2 tests
      fs2.setPartitionAlpha(0.05, false);
      fs2.setSkipEmptyValues(true);

      fs2.findSignificantColumn();

      if (fs2.significantColumnFound()) {
        std::cout << "✓ Best feature: Column " << fs2.getSignificantColumnIndex()
                  << " (" << fs2.getColumnHeader(fs2.getSignificantColumnIndex()) << ")\n";
        std::cout << "  Chi-square: " << fs2.getColumnTestStatistic() << "\n";
        std::cout << "  P-value: " << fs2.getColumnPValue() << "\n";
      } else {
        std::cout << "✗ No significant feature found\n";
      }
    } else {
      std::cout << "(Dataset has fewer than 4 columns; skipping example 2)\n";
    }

    // ========================================================================
    // Example 3: Row filtering (subset of data)
    // ========================================================================
    std::cout << "\n========== Example 3: Row Filtering ==========\n";

    ContingencyTableLib::FeatureSelector fs3;
    fs3.load(outputDir);
    fs3.setTargetColumn(0);

    // Create a row filter: exclude first 10 rows
    std::vector<std::uint32_t> filteredRowsMask(numWords, 0xFFFFFFFF);
    for (std::size_t i = 1; i <= 10 && i < rowCount; ++i) {
      const std::size_t word = i / 32;
      const std::size_t bit = i % 32;
      filteredRowsMask[word] &= ~(1U << bit);
    }

    fs3.enabledRows(filteredRowsMask.data(), rowCount);

    // Enable all candidate columns
    fs3.enabledColumns(candidateColsMask.data(), colCount);

    std::cout << "Testing on rows 11-" << rowCount << "\n";
    std::cout << "Target: Column 0\n";
    std::cout << "Candidates: All except column 0\n\n";

    fs3.setColumnAlpha(0.05, true);
    fs3.setPartitionAlpha(0.05, false);
    fs3.setSkipEmptyValues(true);

    fs3.findSignificantColumn();

    if (fs3.significantColumnFound()) {
      std::cout << "✓ Best feature: Column " << fs3.getSignificantColumnIndex()
                << " (" << fs3.getColumnHeader(fs3.getSignificantColumnIndex()) << ")\n";
      std::cout << "  Chi-square: " << fs3.getColumnTestStatistic() << "\n";
      std::cout << "  P-value: " << fs3.getColumnPValue() << "\n";
    } else {
      std::cout << "✗ No significant feature found\n";
    }

    // ========================================================================
    // Example 4: Stricter alpha with Bonferroni
    // ========================================================================
    std::cout << "\n========== Example 4: Strict Alpha (High Correction) ==========\n";

    ContingencyTableLib::FeatureSelector fs4;
    fs4.load(outputDir);
    fs4.setTargetColumn(0);

    fs4.enabledRows(allRowsMask.data(), rowCount);
    fs4.enabledColumns(candidateColsMask.data(), colCount);

    // Very strict alpha with Bonferroni: 0.001 / (number of columns)
    fs4.setColumnAlpha(0.001, true);  // Stricter threshold
    fs4.setPartitionAlpha(0.05, false);
    fs4.setSkipEmptyValues(true);

    std::cout << "Using strict alpha = 0.001 with Bonferroni correction\n";
    std::cout << "This makes it harder to find significant features\n\n";

    fs4.findSignificantColumn();

    if (fs4.significantColumnFound()) {
      std::cout << "✓ Best feature: Column " << fs4.getSignificantColumnIndex()
                << " (" << fs4.getColumnHeader(fs4.getSignificantColumnIndex()) << ")\n";
      std::cout << "  Chi-square: " << fs4.getColumnTestStatistic() << "\n";
      std::cout << "  P-value: " << fs4.getColumnPValue() << "\n";
    } else {
      std::cout << "✗ No significant feature found at alpha=0.001\n";
    }

    std::cout << "\n========== Summary ==========\n";
    std::cout << "FeatureSelector successfully loaded and analyzed data.\n";
    std::cout << "Key features:\n";
    std::cout << "  • Bonferroni correction for multiple testing\n";
    std::cout << "  • Row filtering via bitmask\n";
    std::cout << "  • Column filtering via bitmask\n";
    std::cout << "  • Configurable alpha levels\n";
    std::cout << "  • Automatic partition search for significant features\n";

    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 2;
  }
}

