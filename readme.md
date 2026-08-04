# MatterGen Scripts

A collection of Python utilities for downloading, analyzing, and preparing crystal structure datasets for MatterGen and other materials science workflows.

## Table of Contents

- [MatterGen Scripts](#mattergen-scripts)
  - [Table of Contents](#table-of-contents)
  - [Features](#features)
  - [Project Structure](#project-structure)
  - [Requirements](#requirements)
  - [Materials Project API Key](#materials-project-api-key)
  - [Scripts](#scripts)
    - [00\_download\_database.py](#00_download_databasepy)
      - [Supported Databases](#supported-databases)
      - [Download Features](#download-features)
      - [Download Example](#download-example)
      - [Download Output](#download-output)
    - [01\_basic\_statistics.py](#01_basic_statisticspy)
      - [Supported File Formats](#supported-file-formats)
      - [Computed Properties](#computed-properties)
      - [Analysis Features](#analysis-features)
      - [Statistics Example](#statistics-example)
      - [Statistics Output](#statistics-output)
  - [Typical Workflow](#typical-workflow)
  - [Output Structure](#output-structure)
  - [Notes](#notes)
  - [License](#license)

---

## Features

- Download crystal structures from supported databases.
- Generate standardized metadata for downloaded datasets.
- Analyze CIF, XYZ, and EXTXYZ crystal structures.
- Compute structural and crystallographic statistics.
- Generate publication-quality figures in PDF and PNG formats.
- Export per-structure and summary CSV files.
- Support parallel processing using all available CPU cores.

---

## Project Structure

```text
.
├── 00_download_database.py
├── 01_basic_statistics.py
├── data/
├── results/
└── README.md
```

---

## Requirements

- Python 3.10+
- pymatgen
- mp-api
- ase
- numpy
- pandas
- scipy
- matplotlib
- seaborn
- tqdm
- joblib
- python-dotenv

Install all dependencies with:

```bash
pip install pymatgen mp-api ase numpy pandas scipy matplotlib seaborn tqdm joblib python-dotenv
```

---

## Materials Project API Key

Downloading from the Materials Project requires an API key.

Create a `.env` file:

```text
MP_API_KEY=YOUR_API_KEY
```

Alternatively, pass the key directly using the `--api-key` argument.

---

## Scripts

### 00_download_database.py

Downloads crystal structures from supported databases and stores them in a standardized directory structure.

#### Supported Databases

| Database          | Status      |
| ----------------- | ----------- |
| Materials Project | Supported   |
| COD               | Placeholder |
| OPTIMADE          | Placeholder |

#### Download Features

- Download one or more chemical systems.
- Save structures as CIF files.
- Automatically create output directories.
- Generate `metadata.csv` for each dataset.
- Skip structures that already exist.
- Print a download summary after completion.

#### Download Example

```bash
python 00_download_database.py \
    --database mp \
    --chemsys Al-O Fe-O Y-Al-O \
    --output ../../data
```

#### Download Output

```text
data/
└── MaterialsProject/
    ├── Al-O/
    │   ├── mp-*.cif
    │   └── metadata.csv
    └── Fe-O/
        ├── mp-*.cif
        └── metadata.csv
```

---

### 01_basic_statistics.py

Computes descriptive statistics for crystal structure datasets.

#### Supported File Formats

- CIF (`.cif`)
- XYZ (`.xyz`)
- EXTXYZ (`.extxyz`)

#### Computed Properties

For every structure, the script computes:

- Chemical formula
- Reduced formula
- Number of atoms
- Number of unique elements
- Lattice parameters (`a`, `b`, `c`)
- Cell angles (`α`, `β`, `γ`)
- Cell volume
- Density
- Space group
- Structure validity

#### Analysis Features

- Parallel processing with Joblib.
- Automatic loading using pymatgen or ASE.
- Summary statistics for all numerical properties.
- Detection of invalid structures.
- CSV export.
- Histogram generation.
- Composition distribution plots.
- Space group distribution plots.
- Lattice parameter distribution plots.

#### Statistics Example

```bash
python 01_basic_statistics.py \
    --input ../../data/MaterialsProject/Al-O \
    --output ../../results/Al-O \
    --workers -1
```

#### Statistics Output

```text
results/
└── Al-O/
    ├── basic_statistics.csv
    ├── summary.csv
    ├── failed_files.csv
    ├── num_atoms_histogram.pdf
    ├── num_atoms_histogram.png
    ├── volume_histogram.pdf
    ├── volume_histogram.png
    ├── density_histogram.pdf
    ├── density_histogram.png
    ├── lattice_parameters.pdf
    ├── lattice_parameters.png
    ├── spacegroup_distribution.pdf
    ├── spacegroup_distribution.png
    ├── composition_distribution.pdf
    └── composition_distribution.png
```

---

## Typical Workflow

1. Download crystal structures from the Materials Project.
2. Store the generated CIF files and metadata.
3. Run the statistics script on the downloaded dataset.
4. Review the generated CSV files.
5. Inspect the generated figures for dataset characteristics.

---

## Output Structure

```text
data/
└── MaterialsProject/
    └── Al-O/
        ├── mp-1000.cif
        ├── mp-1001.cif
        ├── ...
        └── metadata.csv

results/
└── Al-O/
    ├── basic_statistics.csv
    ├── summary.csv
    ├── failed_files.csv
    ├── density_histogram.pdf
    ├── density_histogram.png
    ├── volume_histogram.pdf
    ├── volume_histogram.png
    ├── lattice_parameters.pdf
    ├── lattice_parameters.png
    ├── spacegroup_distribution.pdf
    ├── spacegroup_distribution.png
    ├── composition_distribution.pdf
    └── composition_distribution.png
```

---

## Notes

- Existing downloaded structures are skipped automatically.
- Metadata is generated for every downloaded dataset.
- Invalid structures are listed in `failed_files.csv`.
- Figures are exported as both PDF and PNG.
- The `--workers -1` option uses all available CPU cores.
- Materials Project authentication can be supplied via `.env` or the command line.

---

## License

This project is licensed under the MIT License.
