# AGENTS.md

## Repo purpose (what this library does)
- C++17 library for **contingency tables + Pearson chi-square** on categorical data.
- Provides two user-facing APIs:
  - `ContingencyTableLib::ContingencyTable` for a single pair of columns (+ optional optimal partition search)
  - `ContingencyTableLib::FeatureSelector` for scanning many candidate columns vs a target (Bonferroni-supported) and optionally partitioning the best feature
- Extends / depends on **DataTable v1.1.0** (fetched via CMake `FetchContent`). See `README.md` + `CMakeLists.txt`.

## Big-picture architecture (files to read first)
- Public API headers:
  - `include/ContingencyTable/ContingencyTable.h`
  - `include/ContingencyTable/FeatureSelector.h`
- Implementations:
  - `src/ContingencyTable.cpp`: builds joint counts, computes chi-square + p-value (Wilson–Hilferty), and runs greedy partition search.
  - `src/FeatureSelector.cpp`: orchestrates repeated `ContingencyTable` runs across enabled columns/rows and stores “best” results.
- Examples = canonical usage patterns (keep these working): `examples/*.cpp`.
- Verification notes: `VERIFICATION_RESULTS.md`, `PARTITION_VERIFICATION.md`.

## Data model / integration point (DataTable)
- Inputs are **DataTable “parsed dataset directories”** produced by the `parse_csv` example (`examples/parse_csv.cpp`).
- `ContingencyTable` **inherits** `DataTableLib::DataTable` and reads categorical IDs via `lookupMap(row, col)`.
- Project convention: feature ID **`0` means empty/missing**; when `setSkipEmptyValues(true)` (default) rows with 0 in either column are ignored.

## Core conventions & call ordering (important for correct results)
- Column selection: `setFirstColumn(...)` / `setSecondColumn(...)` accept index or header string (see `ContingencyTable.h`).
- Filtering is via **raw bitmasks** (not owned by the library):
  - `setRowFilter(const uint32_t* bitmask, sizeInBits)` for `ContingencyTable`
  - `enabledRows(...)` / `enabledColumns(...)` for `FeatureSelector`
  - Bit indexes are **0-based**; helpers live inline in `include/ContingencyTable/ContingencyTable.h` (`isBitSet`, `countSetBits`).
- `build()` must be called before reading results; getters throw `std::logic_error` when `dirty_`.
- Partition search: `findOptimalPartition(alpha)` requires a prior `build()` because it uses stored `jointCounts_`.

## Build / test workflows (CMake)
- Configure + build (examples off by default):
  - `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
  - `cmake --build build -j`
- Options in `CMakeLists.txt`:
  - `CONTINGENCYTABLE_BUILD_EXAMPLES` (default OFF)
  - `CONTINGENCYTABLE_BUILD_TESTS` (default ON)
- Run unit tests (CTest): `ctest --test-dir build --output-on-failure`

## Where to add things (project-specific)
- New statistical primitive / table computation: add to `src/ContingencyTable.cpp` and expose in `include/ContingencyTable/ContingencyTable.h`.
- New “scan many columns” behavior: implement in `src/FeatureSelector.cpp`; use `ContingencyTable` as the per-pair engine.
- Add/adjust behavior: update/extend tests in `tests/unit_tests.cpp` (GoogleTest) to lock results; tests also document expected call-order errors.

