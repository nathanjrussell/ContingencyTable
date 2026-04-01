#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "ContingencyTable/ContingencyTable.h"
#include "ContingencyTable/FeatureSelector.h"
#include "DataTable/DataTable.h"

namespace {

std::string makeTempDir(const std::string& suffix) {
  const auto base = std::filesystem::temp_directory_path() /
                    ("contingency_gtest_" + suffix + "_" + std::to_string(std::rand()));
  std::filesystem::create_directories(base);
  return base.string();
}

void cleanupTempDir(const std::string& dir) {
  std::error_code ec;
  std::filesystem::remove_all(dir, ec);
}

std::string writeAndParseCsv(const std::string& tempDir, const std::string& csvContent) {
  const std::string csvPath = tempDir + "/test_data.csv";
  const std::string outputDir = tempDir + "/parsed_output";

  std::ofstream csv(csvPath);
  csv << csvContent;
  csv.close();

  std::filesystem::create_directories(outputDir);

  DataTableLib::DataTable dt(csvPath, outputDir);
  dt.parse(1);

  return outputDir;
}

std::string knownTitanicLikeCsv() {
  return R"(ID,Survived,Pclass
1,No,3
2,No,3
3,No,3
4,No,3
5,No,3
6,No,3
7,No,3
8,No,3
9,No,3
10,No,3
11,No,3
12,No,3
13,No,3
14,No,3
15,No,3
16,No,3
17,No,3
18,No,3
19,No,3
20,No,3
21,No,3
22,No,3
23,No,3
24,No,3
25,No,3
26,No,3
27,No,2
28,No,2
29,No,2
30,No,2
31,No,1
32,No,1
33,No,1
34,No,1
35,Yes,3
36,Yes,3
37,Yes,3
38,Yes,3
39,Yes,3
40,Yes,3
41,Yes,3
42,Yes,3
43,Yes,3
44,Yes,3
45,Yes,2
46,Yes,2
47,Yes,2
48,Yes,2
49,Yes,2
50,Yes,2
51,Yes,2
52,Yes,1
53,Yes,1
54,Yes,1
55,Yes,1
56,Yes,1
57,Yes,1
58,Yes,1
59,Yes,1
60,Yes,1
)";
}

}  // namespace

// Chi-square computation tests
TEST(ContingencyTable, Build_PerfectIndependence_ComputesExpectedChiSquare) {
  const std::string tempDir = makeTempDir("perfect_independence");
  const std::string csv = R"(ID,ColA,ColB
1,X,P
2,X,P
3,X,P
4,X,P
5,Y,Q
6,Y,Q
7,Y,Q
8,Y,Q
)";

  const std::string parseDir = writeAndParseCsv(tempDir, csv);

  ContingencyTableLib::ContingencyTable ct;
  ct.load(parseDir);
  ct.setFirstColumn(1);
  ct.setSecondColumn(2);
  ct.setSkipEmptyValues(true);
  ct.build();

  EXPECT_DOUBLE_EQ(ct.getTestStatistic(), 8.0);
  EXPECT_EQ(ct.getDegreesOfFreedom(), 1);

  cleanupTempDir(tempDir);
}

TEST(ContingencyTable, Build_KnownTable_ComputesExpectedChiSquareAndDf) {
  const std::string tempDir = makeTempDir("known_table");
  const std::string parseDir = writeAndParseCsv(tempDir, knownTitanicLikeCsv());

  ContingencyTableLib::ContingencyTable ct;
  ct.load(parseDir);
  ct.setFirstColumn(1);  // Survived
  ct.setSecondColumn(2); // Pclass
  ct.setSkipEmptyValues(true);
  ct.build();

  EXPECT_NEAR(ct.getTestStatistic(), 8.944720, 1e-3);
  EXPECT_EQ(ct.getDegreesOfFreedom(), 2);

  cleanupTempDir(tempDir);
}

// P-value behavior tests
TEST(ContingencyTable, GetPValue_KnownTable_IsWithinBoundsAndExpectedRange) {
  const std::string tempDir = makeTempDir("pvalue_bounds");
  const std::string parseDir = writeAndParseCsv(tempDir, knownTitanicLikeCsv());

  ContingencyTableLib::ContingencyTable ct;
  ct.load(parseDir);
  ct.setFirstColumn(1);
  ct.setSecondColumn(2);
  ct.setSkipEmptyValues(true);
  ct.build();

  const double p = ct.getPValue();
  EXPECT_GE(p, 0.0);
  EXPECT_LE(p, 1.0);
  EXPECT_NEAR(p, 0.0114, 0.01);

  cleanupTempDir(tempDir);
}

// Partition search tests
TEST(ContingencyTable, FindOptimalPartition_KnownTable_FindsSignificantPartition) {
  const std::string tempDir = makeTempDir("partition_found");
  const std::string parseDir = writeAndParseCsv(tempDir, knownTitanicLikeCsv());

  ContingencyTableLib::ContingencyTable ct;
  ct.load(parseDir);
  ct.setFirstColumn(1);
  ct.setSecondColumn(2);
  ct.setSkipEmptyValues(true);
  ct.build();

  const double originalChi = ct.getTestStatistic();
  auto partition = ct.findOptimalPartition(0.05);

  ASSERT_TRUE(partition.has_value());
  EXPECT_GT(ct.getPartitionChiSquare(), 0.0);
  EXPECT_NE(ct.getPartitionChiSquare(), originalChi);
  EXPECT_LT(ct.getPartitionPValue(), 0.05);
  EXPECT_EQ(ct.getPartitionDegreesOfFreedom(), 1);

  cleanupTempDir(tempDir);
}

// Row filter behavior tests
TEST(ContingencyTable, SetRowFilter_ExcludingRow_ChangesChiSquare) {
  const std::string tempDir = makeTempDir("row_filter");
  const std::string parseDir = writeAndParseCsv(tempDir, knownTitanicLikeCsv());

  ContingencyTableLib::ContingencyTable allRows;
  allRows.load(parseDir);
  allRows.setFirstColumn(1);
  allRows.setSecondColumn(2);
  allRows.setSkipEmptyValues(true);
  allRows.build();
  const double chiAll = allRows.getTestStatistic();

  ContingencyTableLib::ContingencyTable filtered;
  filtered.load(parseDir);
  filtered.setFirstColumn(1);
  filtered.setSecondColumn(2);
  filtered.setSkipEmptyValues(true);

  const auto totalRows = filtered.getRowCount();
  const std::size_t words = (totalRows + 31) / 32;
  std::vector<std::uint32_t> bitmask(words, 0xFFFFFFFF);
  bitmask[0] &= ~(1U << 1);  // Exclude row 1

  filtered.setRowFilter(bitmask.data(), totalRows);
  filtered.build();
  const double chiFiltered = filtered.getTestStatistic();

  EXPECT_NE(chiAll, chiFiltered);
  EXPECT_LT(chiFiltered, chiAll);

  cleanupTempDir(tempDir);
}

// Flag behavior tests
TEST(ContingencyTable, SetSkipEmptyValues_Toggle_ProducesValidStatistics) {
  const std::string tempDir = makeTempDir("skip_empty");
  const std::string csv = R"(ID,ColA,ColB
1,X,P
2,X,Q
3,Y,
4,Y,Q
)";
  const std::string parseDir = writeAndParseCsv(tempDir, csv);

  ContingencyTableLib::ContingencyTable skipEmpty;
  skipEmpty.load(parseDir);
  skipEmpty.setFirstColumn(1);
  skipEmpty.setSecondColumn(2);
  skipEmpty.setSkipEmptyValues(true);
  skipEmpty.build();

  ContingencyTableLib::ContingencyTable includeEmpty;
  includeEmpty.load(parseDir);
  includeEmpty.setFirstColumn(1);
  includeEmpty.setSecondColumn(2);
  includeEmpty.setSkipEmptyValues(false);
  includeEmpty.build();

  EXPECT_GE(skipEmpty.getTestStatistic(), 0.0);
  EXPECT_GE(includeEmpty.getTestStatistic(), 0.0);

  cleanupTempDir(tempDir);
}

TEST(FeatureSelector, GetTargetCounts_CachesMarginalsAndClearsOnRerun) {
  const std::string tempDir = makeTempDir("feature_selector_target_counts");
  // Target column (T) has A,A,A,B,B across included rows.
  // Candidate column (X) matches T perfectly so it should be selected as best.
  const std::string csv = R"(ID,T,X
1,A,A
2,A,A
3,A,A
4,B,B
5,B,B
)";

  const std::string parseDir = writeAndParseCsv(tempDir, csv);

  ContingencyTableLib::FeatureSelector fs;
  fs.load(parseDir);
  fs.setTargetColumn(1);

  // Enable all rows.
  const auto rowCount = static_cast<std::size_t>(fs.getRowCount());
  const std::size_t rowWords = (rowCount + 31) / 32;
  std::vector<std::uint32_t> allRowsMask(rowWords, 0xFFFFFFFF);
  fs.enabledRows(allRowsMask.data(), rowCount);

  // Enable only column 2 as candidate (X).
  const auto colCount = static_cast<std::size_t>(fs.getColumnCount());
  std::vector<std::uint32_t> colMask((colCount + 31) / 32, 0);
  colMask[0] |= (1U << 2);
  fs.enabledColumns(colMask.data(), colCount);

  fs.setSkipEmptyValues(true);
  fs.setColumnAlpha(0.5, false);  // Make it easy to be "significant" in small test data.
  fs.setPartitionAlpha(0.5, false);

  fs.findSignificantColumn();
  ASSERT_TRUE(fs.significantColumnFound());
  EXPECT_EQ(fs.getSignificantColumnIndex(), 2U);

  // Compute expected mapped IDs from DataTable.
  DataTableLib::DataTable dt;
  dt.load(parseDir);
  const auto idA = dt.lookupMap(1, 1);
  const auto idB = dt.lookupMap(4, 1);

  const auto counts = fs.getTargetCounts();
  ASSERT_TRUE(counts.find(idA) != counts.end());
  ASSERT_TRUE(counts.find(idB) != counts.end());
  EXPECT_EQ(counts.at(idA), 3ULL);
  EXPECT_EQ(counts.at(idB), 2ULL);

  // Ensure cache is wiped on rerun-triggering configuration changes.
  fs.enabledRows(allRowsMask.data(), rowCount);  // resets and clears cached data
  fs.enabledColumns(colMask.data(), colCount);
  fs.findSignificantColumn();
  ASSERT_TRUE(fs.significantColumnFound());
  const auto counts2 = fs.getTargetCounts();
  EXPECT_EQ(counts2.at(idA), 3ULL);
  EXPECT_EQ(counts2.at(idB), 2ULL);

  cleanupTempDir(tempDir);
}

// ============================================================================
// FeatureSelector Tests
// ============================================================================

TEST(FeatureSelector, FindSignificantColumn_WithBonferroniCorrection_FindsStrongSignal) {
  const std::string tempDir = makeTempDir("fs_bonferroni");

  // Create dataset with target and 3 candidate features
  // feature2 has strong association with target
  std::string csv = "target,feature1,feature2,feature3\n";
  for (int i = 0; i < 100; ++i) {
    csv += (i < 50 ? "A" : "B");
    csv += ",";
    csv += (i % 2 == 0 ? "X" : "Y");  // weak
    csv += ",";
    csv += (i < 50 ? "P" : "Q");      // strong
    csv += ",";
    csv += (i % 3 == 0 ? "M" : "N");  // weak
    csv += "\n";
  }

  const std::string parseDir = writeAndParseCsv(tempDir, csv);

  ContingencyTableLib::FeatureSelector fs;
  fs.load(parseDir);
  fs.setTargetColumn(0);

  // Enable all rows and columns 1, 2, 3
  std::vector<std::uint32_t> rowMask(4, 0xFFFFFFFF);
  std::uint32_t colMask = 0b1110; // bits 1, 2, 3 enabled
  fs.enabledRows(rowMask.data(), 128);
  fs.enabledColumns(&colMask, 4);

  // With Bonferroni correction
  fs.setColumnAlpha(0.05, true);
  fs.findSignificantColumn();

  EXPECT_TRUE(fs.significantColumnFound());
  EXPECT_EQ(fs.getSignificantColumnIndex(), 2);
  EXPECT_GT(fs.getColumnTestStatistic(), 0.0);
  EXPECT_LT(fs.getColumnPValue(), 0.05 / 3); // Bonferroni adjusted

  cleanupTempDir(tempDir);
}

TEST(FeatureSelector, ThrowsWhenAccessingResultsBeforeSearch) {
  ContingencyTableLib::FeatureSelector fs;

  EXPECT_THROW(fs.getSignificantColumnIndex(), std::runtime_error);
  EXPECT_THROW(fs.getColumnPValue(), std::runtime_error);
  EXPECT_THROW(fs.getColumnTestStatistic(), std::runtime_error);
  EXPECT_THROW(fs.getColumnDegreesOfFreedom(), std::runtime_error);
  EXPECT_THROW(fs.getFirstPartition(), std::runtime_error);
  EXPECT_THROW(fs.getSecondPartition(), std::runtime_error);
}

TEST(FeatureSelector, NoSignificantColumnFound_ReturnsFalse) {
  const std::string tempDir = makeTempDir("fs_no_signal");

  // Random noise - no correlation
  std::string csv = "target,noise\n";
  for (int i = 0; i < 50; ++i) {
    csv += (i % 2 == 0 ? "A" : "B");
    csv += ",";
    csv += (i % 3 == 0 ? "X" : "Y");
    csv += "\n";
  }

  const std::string parseDir = writeAndParseCsv(tempDir, csv);

  ContingencyTableLib::FeatureSelector fs;
  fs.load(parseDir);
  fs.setTargetColumn(0);

  std::vector<std::uint32_t> rowMask(2, 0xFFFFFFFF);
  std::uint32_t colMask = 0b10; // only column 1
  fs.enabledRows(rowMask.data(), 64);
  fs.enabledColumns(&colMask, 2);

  fs.setColumnAlpha(0.001, false);  // Very strict threshold
  fs.findSignificantColumn();

  EXPECT_FALSE(fs.significantColumnFound());
  EXPECT_FALSE(fs.significantPartitionFound());

  cleanupTempDir(tempDir);
}

TEST(FeatureSelector, PartitionSearchOnSignificantColumn_FindsPartition) {
  const std::string tempDir = makeTempDir("fs_partition");
  const std::string parseDir = writeAndParseCsv(tempDir, knownTitanicLikeCsv());

  ContingencyTableLib::FeatureSelector fs;
  fs.load(parseDir);
  fs.setTargetColumn(1);  // Survived

  std::vector<std::uint32_t> rowMask(2, 0xFFFFFFFF);
  std::uint32_t colMask = 0b100; // only column 2 (Pclass)
  fs.enabledRows(rowMask.data(), 64);
  fs.enabledColumns(&colMask, 3);

  fs.setPartitionAlpha(0.05, false);
  fs.findSignificantColumn();

  ASSERT_TRUE(fs.significantColumnFound());
  EXPECT_EQ(fs.getSignificantColumnIndex(), 2);

  if (fs.significantPartitionFound()) {
    auto p0 = fs.getFirstPartition();
    auto p1 = fs.getSecondPartition();

    EXPECT_FALSE(p0.empty());
    EXPECT_FALSE(p1.empty());
    EXPECT_GT(fs.getPartitionTestStatistic(), 0.0);
    EXPECT_LT(fs.getPartitionPValue(), 0.05);
    EXPECT_EQ(fs.getPartitionDegreesOfFreedom(), 1);
  }

  cleanupTempDir(tempDir);
}

TEST(FeatureSelector, PartitionRowCounts_CountIncludedRowsAndRespectRowFilter) {
  const std::string tempDir = makeTempDir("fs_partition_row_counts");

  // Reuse the known Titanic-like dataset which is already used elsewhere in this suite.
  // It reliably produces a significant association between Survived (col 1) and Pclass (col 2).
  const std::string parseDir = writeAndParseCsv(tempDir, knownTitanicLikeCsv());

  ContingencyTableLib::FeatureSelector fs;
  fs.load(parseDir);
  fs.setTargetColumn(1); // Survived
  fs.setSkipEmptyValues(true);

  // Only test Pclass as candidate so numTests=1 and significance matches the CT-based checks.
  std::uint32_t colMask = 0b100; // enable column 2 only
  fs.enabledColumns(&colMask, 3);

  const auto rowCount = fs.getRowCount();
  const std::size_t rowWords = (static_cast<std::size_t>(rowCount) + 31) / 32;

  // Helper to compute included row count for a given row mask.
  auto computeIncludedRows = [&](const std::vector<std::uint32_t>& rowMask) -> std::uint64_t {
    ContingencyTableLib::ContingencyTable ct;
    ct.load(parseDir);
    ct.setFirstColumn(1);
    ct.setSecondColumn(2);
    ct.setSkipEmptyValues(true);
    ct.setRowFilter(rowMask.data(), static_cast<std::size_t>(rowCount));
    ct.build();

    // Included rows are the sum over all jointCounts_, which we can obtain by partitioning
    // all observed B values into partition0 (all of them) vs partition1 (none).
    std::vector<std::uint32_t> allB;
    for (std::uint64_t r = 1; r < ct.getRowCount(); ++r) {
      const auto b = ct.lookupMap(r, 2);
      const auto a = ct.lookupMap(r, 1);
      if (a == 0 || b == 0) {
        continue;
      }
      allB.push_back(b);
    }
    // allB may contain duplicates; de-dup to keep the helper cheap.
    std::sort(allB.begin(), allB.end());
    allB.erase(std::unique(allB.begin(), allB.end()), allB.end());
    const auto counts = ct.countRowsInPartitions(allB);
    return counts.first + counts.second;
  };

  // Enable all rows initially.
  std::vector<std::uint32_t> allRowsMask(rowWords, 0xFFFFFFFF);
  fs.enabledRows(allRowsMask.data(), static_cast<std::size_t>(rowCount));

  fs.findSignificantColumn();

  ASSERT_TRUE(fs.significantColumnFound());
  EXPECT_EQ(fs.getSignificantColumnIndex(), 2U);

  if (fs.significantPartitionFound()) {
    const auto total = fs.getPartitionOneRowCount() + fs.getPartitionTwoRowCount();
    EXPECT_EQ(total, computeIncludedRows(allRowsMask));

    // Now exclude one data row via row filter and ensure the sum decreases by 1.
    // Exclude data row index 1 (first data row) in DataTable indexing.
    std::vector<std::uint32_t> filteredMask(rowWords, 0xFFFFFFFF);
    filteredMask[0] &= ~(1U << 1);
    fs.enabledRows(filteredMask.data(), static_cast<std::size_t>(rowCount));
    fs.findSignificantColumn();

    ASSERT_TRUE(fs.significantColumnFound());
    if (fs.significantPartitionFound()) {
      const auto filteredTotal = fs.getPartitionOneRowCount() + fs.getPartitionTwoRowCount();
      EXPECT_EQ(filteredTotal, computeIncludedRows(filteredMask));
    }
  }

  cleanupTempDir(tempDir);
}

TEST(FeatureSelector, SkipEmptyValues_IsRespected) {
  const std::string tempDir = makeTempDir("fs_empty");

  std::string csv = "target,feature\n";
  csv += "A,X\n";
  csv += "A,\n";  // empty value
  csv += "B,Y\n";
  csv += "B,X\n";

  const std::string parseDir = writeAndParseCsv(tempDir, csv);

  ContingencyTableLib::FeatureSelector fs;
  fs.load(parseDir);
  fs.setTargetColumn(0);
  fs.setSkipEmptyValues(true);

  std::vector<std::uint32_t> rowMask(1, 0xFFFFFFFF);
  std::uint32_t colMask = 0b10;
  fs.enabledRows(rowMask.data(), 32);
  fs.enabledColumns(&colMask, 2);

  EXPECT_NO_THROW(fs.findSignificantColumn());

  cleanupTempDir(tempDir);
}

TEST(FeatureSelector, BonferroniCorrection_AdjustsAlphaCorrectly) {
  const std::string tempDir = makeTempDir("fs_bonf_adjust");

  // Create dataset with many weak features and one strong feature
  std::string csv = "target,f1,f2,f3,f4,f5,strong\n";
  for (int i = 0; i < 100; ++i) {
    csv += (i < 50 ? "A" : "B");
    for (int j = 0; j < 5; ++j) {
      csv += "," + std::string(1, 'X' + (i % 2));  // weak noise
    }
    csv += "," + std::string(i < 50 ? "P" : "Q");  // strong signal
    csv += "\n";
  }

  const std::string parseDir = writeAndParseCsv(tempDir, csv);

  ContingencyTableLib::FeatureSelector fs;
  fs.load(parseDir);
  fs.setTargetColumn(0);

  std::vector<std::uint32_t> rowMask(4, 0xFFFFFFFF);
  std::uint32_t colMask = 0b1111110; // columns 1-6
  fs.enabledRows(rowMask.data(), 128);
  fs.enabledColumns(&colMask, 7);

  fs.setColumnAlpha(0.05, true);  // Bonferroni with 6 tests
  fs.findSignificantColumn();

  // Should find the strong feature even with Bonferroni correction
  EXPECT_TRUE(fs.significantColumnFound());
  EXPECT_EQ(fs.getSignificantColumnIndex(), 6);

  cleanupTempDir(tempDir);
}

TEST(FeatureSelector, RowFiltering_AffectsResults) {
  const std::string tempDir = makeTempDir("fs_row_filter");
  const std::string parseDir = writeAndParseCsv(tempDir, knownTitanicLikeCsv());

  // Test with all rows
  ContingencyTableLib::FeatureSelector fs1;
  fs1.load(parseDir);
  fs1.setTargetColumn(1);

  std::vector<std::uint32_t> allRowsMask(2, 0xFFFFFFFF);
  std::uint32_t colMask = 0b100; // column 2
  fs1.enabledRows(allRowsMask.data(), 64);
  fs1.enabledColumns(&colMask, 3);
  fs1.findSignificantColumn();

  const double pval1 = fs1.significantColumnFound() ? fs1.getColumnPValue() : 1.0;

  // Test with filtered rows
  ContingencyTableLib::FeatureSelector fs2;
  fs2.load(parseDir);
  fs2.setTargetColumn(1);

  std::vector<std::uint32_t> filteredMask(2, 0xFFFFFFFF);
  filteredMask[0] &= ~(1U << 1);  // exclude row 1
  filteredMask[0] &= ~(1U << 2);  // exclude row 2
  fs2.enabledRows(filteredMask.data(), 64);
  fs2.enabledColumns(&colMask, 3);
  fs2.findSignificantColumn();

  const double pval2 = fs2.significantColumnFound() ? fs2.getColumnPValue() : 1.0;

  // P-values should differ when rows are filtered
  EXPECT_NE(pval1, pval2);

  cleanupTempDir(tempDir);
}
