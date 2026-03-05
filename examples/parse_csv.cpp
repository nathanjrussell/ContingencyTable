#include <DataTable/DataTable.h>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Usage: parse_csv <input.csv> <output_dir>\n";
    return 1;
  }

  try {
    const std::string inputCsv = argv[1];
    const std::string outputDir = argv[2];

    DataTableLib::DataTable table(inputCsv, outputDir);

    std::cout << "Parsing CSV: " << inputCsv << "\n";
    table.parse(4);  // Use 4 threads

    std::cout << "Successfully parsed!\n";
    std::cout << "Output directory: " << outputDir << "\n";
    std::cout << "Rows: " << table.getRowCount() << "\n";
    std::cout << "Columns: " << table.getColumnCount() << "\n";

    std::cout << "\nColumn headers:\n";
    for (std::uint64_t i = 0; i < table.getColumnCount(); ++i) {
      std::cout << "  " << i << ": " << table.getColumnHeader(i) << "\n";
    }

    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 2;
  }
}

