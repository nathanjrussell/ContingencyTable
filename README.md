# ContingencyTable

A C++ library for computing contingency tables and chi-square statistics from categorical data. Extends [DataTable](https://github.com/nathanjrussell/DataTable) v1.1.0.

## Features

- **ContingencyTable**: Column selection by index or header name, row filtering via bitmask, chi-square test with p-value (Wilson-Hilferty approximation), optimal partition search
- **FeatureSelector**: Automated feature selection across multiple columns with Bonferroni correction, finds most significant feature and optimal partition
- No external dependencies beyond DataTable
- SLURM-friendly (no Boost required)

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCONTINGENCYTABLE_BUILD_EXAMPLES=ON
cmake --build build -j
```

## Usage

### 1. Parse CSV

```bash
./build/parse_csv  output_dir
```

### 2. Run Chi-Square Analysis

```bash
./build/contingency_example output_dir col1_index col2_index
```

## Example: Titanic Dataset

```bash
# Parse the data
./build/parse_csv examples/datasets/data_set_1.csv examples/datasets/titanic_output

# Test: Survived vs Sex
./build/contingency_example examples/datasets/titanic_output 1 3
```

Output:
```
Loaded DataTable: 61 rows, 8 columns
Chi-square statistic: 15.0179
Degrees of freedom:   1
P-value:              0.000169323
```

Highly significant (p < 0.001) - females had much higher survival rates.

```bash
# Test: Survived vs Passenger Class
./build/contingency_example examples/datasets/titanic_output 1 2
```

Output:
```
Loaded DataTable: 61 rows, 8 columns
Chi-square statistic: 8.94472
Degrees of freedom:   2
P-value:              0.0114211
```

Significant (p < 0.05) - first-class passengers had better survival rates.

### 3. Run Feature Selection

```bash
./build/feature_selector examples/datasets/titanic_output
```

This demonstrates automated feature selection with the FeatureSelector class:

Output:
```
Loaded DataTable from: examples/datasets/titanic_output
Total rows: 61
Total columns: 8

========== Example 1: Feature Selection ==========
Target: Column 0 (first column)
Candidates: All other columns

Searching for most significant feature...
✗ No significant feature found at alpha=0.05

========== Example 2: Select Specific Columns ==========
Target: Column 1 (Survived)
Candidates: Column 2 (Pclass), Column 3 (Sex)

✓ Best feature: Column 3 (Sex)
  Chi-square: 15.0179
  P-value: 0.000169323

========== Example 3: Row Filtering ==========
Testing on rows 11-61
Target: Column 0
Candidates: All except column 0

✗ No significant feature found

========== Example 4: Strict Alpha (High Correction) ==========
Using strict alpha = 0.001 with Bonferroni correction
This makes it harder to find significant features

✗ No significant feature found at alpha=0.001
```

The FeatureSelector automatically:
- Tests all candidate columns against the target
- Applies Bonferroni correction for multiple testing
- Identifies the most significant feature (Sex in Example 2)
- Can search for optimal partitions on significant features

## API

### ContingencyTable

```cpp
#include <ContingencyTable/ContingencyTable.h>
#include <iostream>
#include <vector>

int main() {
    ContingencyTableLib::ContingencyTable table;
    
    // Load the parsed Titanic dataset
    table.load("examples/datasets/titanic_output");
    
    // Select columns by index or name
    table.setFirstColumn(1);        // Survived column
    table.setSecondColumn("Sex");   // Sex column (by name)
    
    // Optional: row filtering using bitmask
    const auto rowCount = table.getRowCount();
    const std::size_t numWords = (rowCount + 31) / 32;
    std::vector<std::uint32_t> bitmask(numWords, 0xFFFFFFFF);
    
    // Example: skip row 5
    bitmask[5 / 32] &= ~(1U << (5 % 32));
    table.setRowFilter(bitmask.data(), rowCount);
    
    // Optional: empty value handling
    table.setSkipEmptyValues(true);  // Default: skip rows with 0 values
    
    // Compute the contingency table and statistics
    table.build();
    
    // Get and display results
    std::cout << "Chi-square: " << table.getTestStatistic() << "\n";
    std::cout << "DF: " << table.getDegreesOfFreedom() << "\n";
    std::cout << "P-value: " << table.getPValue() << "\n";
    
    return 0;
}
```

### FeatureSelector

```cpp
#include <ContingencyTable/FeatureSelector.h>
#include <iostream>
#include <vector>

int main() {
    ContingencyTableLib::FeatureSelector fs;
    
    // Load parsed dataset
    fs.load("examples/datasets/titanic_output");
    fs.setTargetColumn(1);  // Survived column
    
    // Enable rows (all rows in this example)
    std::vector<std::uint32_t> rowMask(4, 0xFFFFFFFF);
    fs.enabledRows(rowMask.data(), 128);
    
    // Enable candidate columns (columns 2, 3, 4 via bitmask)
    std::uint32_t colMask = 0b11100;  // bits 2, 3, 4 set
    fs.enabledColumns(&colMask, 5);
    
    // Configure alpha thresholds with Bonferroni correction
    fs.setColumnAlpha(0.05, true);      // Apply Bonferroni for column selection
    fs.setPartitionAlpha(0.05, false);  // No correction for partition search
    fs.setSkipEmptyValues(true);
    
    // Find most significant column
    fs.findSignificantColumn();
    
    if (fs.significantColumnFound()) {
        std::cout << "Best feature: column " << fs.getSignificantColumnIndex() << "\n";
        std::cout << "Chi-square: " << fs.getColumnTestStatistic() << "\n";
        std::cout << "P-value: " << fs.getColumnPValue() << "\n";
        std::cout << "DF: " << fs.getColumnDegreesOfFreedom() << "\n";
        
        // Check if a partition was found
        if (fs.significantPartitionFound()) {
            auto p0 = fs.getFirstPartition();
            auto p1 = fs.getSecondPartition();
            std::cout << "Partition found with " << p0.size() 
                      << " and " << p1.size() << " categories\n";
            std::cout << "Partition chi-square: " << fs.getPartitionTestStatistic() << "\n";
        }
    } else {
        std::cout << "No significant feature found\n";
    }
    
    return 0;
}
```

**To compile and run this example:**
```bash
# First, parse the Titanic dataset
./build/parse_csv examples/datasets/data_set_1.csv examples/datasets/titanic_output

# Compile your program
g++ -std=c++17 -I./include -L./build -o my_analysis my_analysis.cpp \
    -lContingencyTable -lDataTable

# Run it
./my_analysis
```

## Validation

Chi-square values match scipy.stats.chi2_contingency with `correction=False`. Row filtering verified by comparing filtered results between C++ and scipy.

**Note:** scipy uses Yates' correction by default for 2×2 tables, which reduces chi-square values. Use `correction=False` in scipy to match this library.

## Project Structure

```
include/ContingencyTable/ContingencyTable.h  - Public API
src/ContingencyTable.cpp                      - Implementation
examples/parse_csv.cpp                        - CSV parser
examples/contingency_example.cpp              - Chi-square example
examples/datasets/data_set_1.csv              - Sample Titanic data
```

## Dependencies

- DataTable v1.1.0 (fetched via CMake FetchContent)
- C++17 compiler
- CMake 3.21+

## Unit Tests

Unit tests use Google Test and are registered with CTest.

### Build With Tests

Tests are enabled by default (`CONTINGENCYTABLE_BUILD_TESTS=ON`):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Run All Tests

```bash
ctest --test-dir build --output-on-failure
```

### Run Test Binary Directly (Optional)

```bash
./build/ContingencyTableTests
./build/ContingencyTableTests --gtest_list_tests
./build/ContingencyTableTests --gtest_filter=FeatureSelector.*
```

### Current Test Cases (13 total)

**ContingencyTable (6 tests):**
- Build_PerfectIndependence_ComputesExpectedChiSquare
- Build_KnownTable_ComputesExpectedChiSquareAndDf
- GetPValue_KnownTable_IsWithinBoundsAndExpectedRange
- FindOptimalPartition_KnownTable_FindsSignificantPartition
- SetRowFilter_ExcludingRow_ChangesChiSquare
- SetSkipEmptyValues_Toggle_ProducesValidStatistics

**FeatureSelector (7 tests):**
- FindSignificantColumn_WithBonferroniCorrection_FindsStrongSignal
- ThrowsWhenAccessingResultsBeforeSearch
- NoSignificantColumnFound_ReturnsFalse
- PartitionSearchOnSignificantColumn_FindsPartition
- SkipEmptyValues_IsRespected
- BonferroniCorrection_AdjustsAlphaCorrectly
- RowFiltering_AffectsResults

### Notes

- Tests create temporary CSV datasets under `/tmp` and clean up automatically.
- To build examples and tests together:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCONTINGENCYTABLE_BUILD_EXAMPLES=ON -DCONTINGENCYTABLE_BUILD_TESTS=ON
cmake --build build -j
```
