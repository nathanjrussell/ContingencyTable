# ContingencyTable

A C++ library for computing contingency tables and chi-square statistics from categorical data. Extends [DataTable](https://github.com/nathanjrussell/DataTable) v1.1.0.

## Features

- Column selection by index or header name
- Row filtering via bitmask
- Chi-square test with p-value (Wilson-Hilferty approximation)
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

## API

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



