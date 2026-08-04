# MatterGen Scripts

A collection of Python scripts for downloading, analyzing, and benchmarking crystal structure datasets from MatterGen and established materials databases.

The goal of this project is to evaluate the quality, realism, diversity, and novelty of MatterGen-generated crystal structures by comparing them against large experimental and computational crystal databases using a unified analysis pipeline.

---

## Table of Contents

- [MatterGen Scripts](#mattergen-scripts)
  - [Table of Contents](#table-of-contents)
  - [Features](#features)
  - [Project Structure](#project-structure)
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Materials Project API Key](#materials-project-api-key)
  - [Scripts](#scripts)
    - [00\_download\_database.py](#00_download_databasepy)
      - [Download Features](#download-features)
      - [Download Example](#download-example)
    - [01\_basic\_statistics.py](#01_basic_statisticspy)
      - [Supported File Formats](#supported-file-formats)
      - [Computed Properties](#computed-properties)
      - [Statistics Output](#statistics-output)
      - [Statistics Example](#statistics-example)
    - [02\_symmetry\_analysis.py](#02_symmetry_analysispy)
      - [Symmetry Features](#symmetry-features)
      - [Symmetry: Computed Properties](#symmetry-computed-properties)
      - [Symmetry Output](#symmetry-output)
      - [Symmetry Example](#symmetry-example)
    - [03\_coordination\_analysis.py](#03_coordination_analysispy)
      - [Coordination Features](#coordination-features)
      - [Coordination: Computed Properties](#coordination-computed-properties)
      - [Coordination Output](#coordination-output)
      - [Coordination Example](#coordination-example)
  - [Implemented Analysis Pipeline](#implemented-analysis-pipeline)
  - [Planned Analysis Pipeline](#planned-analysis-pipeline)
  - [Directory Structure](#directory-structure)
  - [Typical Workflow](#typical-workflow)
  - [Supported Databases](#supported-databases)
  - [Roadmap](#roadmap)
  - [License](#license)

---

## Features

- Download crystal structures from multiple materials databases.
- Support Materials Project and OPTIMADE-compatible databases.
- Download multiple chemical systems in a single run.
- Store all downloaded structures as CIF files.
- Automatically create the required directory structure.
- Skip already downloaded structures.
- Analyze CIF, XYZ, and EXTXYZ crystal structures.
- Compute structural and crystallographic statistics.
- Analyze crystallographic symmetry.
- Determine space groups, crystal systems, point groups, Bravais lattices, and Hall symbols.
- Analyze local coordination environments using CrystalNN and VoronoiNN.
- Export CSV summaries.
- Generate publication-quality PDF and PNG figures.
- Support parallel processing using all available CPU cores.
- Provide a common benchmarking workflow for MatterGen and reference datasets.

---

## Project Structure

```text
.
├── 00_download_database.py
├── 01_basic_statistics.py
├── 02_symmetry_analysis.py
├── 03_coordination_analysis.py
├── downloaders/
│   ├── __init__.py
│   ├── alexandria.py
│   ├── base.py
│   ├── cod.py
│   ├── jarvis.py
│   ├── materials_project.py
│   ├── oqmd.py
│   ├── optimade.py
│   └── utils.py
├── data/
├── results/
├── .env.example
├── .gitignore
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
- requests
- python-dotenv

---

## Installation

Clone the repository:

```bash
git clone https://github.com/<username>/mattergen-scripts.git

cd mattergen-scripts
```

Install the required packages:

```bash
pip install \
    pymatgen \
    mp-api \
    ase \
    numpy \
    pandas \
    scipy \
    matplotlib \
    seaborn \
    tqdm \
    joblib \
    requests \
    python-dotenv
```

---

## Materials Project API Key

Materials Project downloads require an API key.

Create a `.env` file by copying the example:

Linux/macOS

```bash
cp .env.example .env
```

Windows PowerShell

```powershell
Copy-Item .env.example .env
```

Edit the `.env` file:

```text
MP_API_KEY=YOUR_API_KEY
```

The downloader automatically loads the API key from the `.env` file.

You can also override it manually:

```bash
python 00_download_database.py \
    --database mp \
    --api-key YOUR_API_KEY
```

---

## Scripts

### 00_download_database.py

Downloads crystal structures from supported databases and stores them as CIF files.

#### Download Features

- Download one or more chemical systems
- Download multiple databases
- Save structures as CIF files
- Skip existing downloads
- Automatically create output folders
- Progress bars
- Download summaries

#### Download Example

```bash
python 00_download_database.py \
    --database mp \
    --chemsys Al-O Fe-O Ti-O Y-Al-O \
    --output ../../data
```

---

### 01_basic_statistics.py

Computes descriptive statistics for crystal structure datasets.

#### Supported File Formats

- CIF (`.cif`)
- XYZ (`.xyz`)
- EXTXYZ (`.extxyz`)

#### Computed Properties

For every structure:

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

#### Statistics Output

Results are written to a dedicated `basic_statistics` directory.

```text
results/
└── <database>/
    └── <chemical-system>/
        └── basic_statistics/
            ├── basic_statistics.csv
            ├── summary.csv
            ├── failed_files.csv
            ├── composition_distribution.pdf
            ├── composition_distribution.png
            ├── density_histogram.pdf
            ├── density_histogram.png
            ├── lattice_parameters.pdf
            ├── lattice_parameters.png
            ├── num_atoms_histogram.pdf
            ├── num_atoms_histogram.png
            ├── spacegroup_distribution.pdf
            ├── spacegroup_distribution.png
            ├── volume_histogram.pdf
            └── volume_histogram.png
```

#### Statistics Example

```bash
python 01_basic_statistics.py \
    --input ../../data/MaterialsProject/Al-O \
    --output ../../results/MaterialsProject/Al-O \
    --workers -1
```

---

### 02_symmetry_analysis.py

Analyzes the crystallographic symmetry of crystal structure datasets.

#### Symmetry Features

- Analyze CIF, XYZ, and EXTXYZ crystal structures
- Determine crystallographic symmetry using `SpacegroupAnalyzer`
- Compute one symmetry analysis per structure
- Export CSV summaries
- Generate publication-quality PDF and PNG figures
- Support parallel processing

#### Symmetry: Computed Properties

For every structure:

- Space group number
- Space group symbol
- Crystal system
- Bravais lattice
- Point group
- Hall symbol
- Number of symmetry operations
- Structure validity

#### Symmetry Output

Results are written to a dedicated `symmetry` directory.

```text
results/
└── <database>/
    └── <chemical-system>/
        └── symmetry/
            ├── symmetry_analysis.csv
            ├── summary.csv
            ├── failed_files.csv
            ├── crystal_system_distribution.pdf
            ├── crystal_system_distribution.png
            ├── point_group_distribution.pdf
            ├── point_group_distribution.png
            ├── space_group_distribution.pdf
            ├── space_group_distribution.png
            ├── space_group_number_histogram.pdf
            ├── space_group_number_histogram.png
            ├── symmetry_operations_histogram.pdf
            └── symmetry_operations_histogram.png
```

#### Symmetry Example

```bash
python 02_symmetry_analysis.py \
    --input ../../data/MaterialsProject/Al-O \
    --output ../../results/MaterialsProject/Al-O \
    --workers -1
```

---

### 03_coordination_analysis.py

Analyzes local coordination environments of crystal structure datasets using both CrystalNN and VoronoiNN.

#### Coordination Features

- Analyze CIF, XYZ, and EXTXYZ crystal structures
- Compute local coordination using CrystalNN
- Compute local coordination using VoronoiNN
- Compare coordination assignments from both methods
- Compute per-structure coordination statistics
- Export CSV summaries
- Generate publication-quality PDF and PNG figures
- Support parallel processing

#### Coordination: Computed Properties

For every structure:

- Mean CrystalNN coordination number
- Median CrystalNN coordination number
- Standard deviation of CrystalNN coordination
- Minimum CrystalNN coordination
- Maximum CrystalNN coordination
- Mean VoronoiNN coordination number
- Median VoronoiNN coordination number
- Standard deviation of VoronoiNN coordination
- Minimum VoronoiNN coordination
- Maximum VoronoiNN coordination
- Mean absolute coordination difference
- Maximum absolute coordination difference
- Agreement fraction between both methods
- Structure validity

#### Coordination Output

Results are written to a dedicated `coordination` directory.

```text
results/
└── <database>/
    └── <chemical-system>/
        └── coordination/
            ├── coordination_analysis.csv
            ├── summary.csv
            ├── failed_files.csv
            ├── agreement_fraction_histogram.pdf
            ├── agreement_fraction_histogram.png
            ├── composition_distribution.pdf
            ├── composition_distribution.png
            ├── coordination_boxplot.pdf
            ├── coordination_boxplot.png
            ├── coordination_difference_histogram.pdf
            ├── coordination_difference_histogram.png
            ├── coordination_scatter.pdf
            ├── coordination_scatter.png
            ├── crystalnn_coordination_histogram.pdf
            ├── crystalnn_coordination_histogram.png
            ├── voronoi_coordination_histogram.pdf
            └── voronoi_coordination_histogram.png
```

#### Coordination Example

```bash
python 03_coordination_analysis.py     --input ../../data/MaterialsProject/Al-O     --output ../../results/MaterialsProject/Al-O     --workers -1
```

## Implemented Analysis Pipeline

| Script | Analysis |
| ------ | -------- |
| 00_download_database.py | Download reference datasets |
| 01_basic_statistics.py | Basic structural statistics |
| 02_symmetry_analysis.py | Space groups, crystal systems, Bravais lattices, point groups, Hall symbols, and symmetry operations |
| 03_coordination_analysis.py | CrystalNN and VoronoiNN coordination environments, agreement metrics, and coordination statistics |

---

## Planned Analysis Pipeline

| Script | Analysis |
| ------ | -------- |
| 04_bond_analysis.py | Bond lengths, bond angles, nearest neighbors |
| 05_rdf_analysis.py | Radial distribution functions (RDF) |
| 06_stability_analysis.py | Formation energy, energy above hull, band gap and stability metrics |

---

## Directory Structure

```text
data/
├── MatterGen/
│   ├── Al-O/
│   │   └── generated_crystals_cif/
│   ├── Fe-O/
│   ├── Ti-O/
│   └── Y-Al-O/
├── MaterialsProject/
├── COD/
├── OQMD/
├── JARVIS/
└── Alexandria/

results/
├── MatterGen/
│   └── Al-O/
│       ├── basic_statistics/
│       ├── symmetry/
│       └── coordination/
├── MaterialsProject/
├── COD/
├── OQMD/
├── JARVIS/
└── Alexandria/
```

---

## Typical Workflow

1. Download reference crystal structures from one or more databases.
2. Generate crystal structures with MatterGen.
3. Run the implemented analysis pipeline on every dataset.
4. Compare structural statistics.
5. Compare crystallographic symmetry.
6. Compare local coordination environments.
7. Continue with the remaining planned analyses.
8. Determine whether MatterGen generates realistic crystal structures.

---

## Supported Databases

| Database | Type | Status |
| -------- | ------------ | :----: |
| Materials Project | DFT | ✅ |
| COD | Experimental | ✅ |
| JARVIS | DFT | ✅ |
| OQMD | DFT | ✅ |
| Alexandria | DFT / ML | ⚠ Experimental |

---

## Roadmap

The benchmark aims to compare:

- Chemical composition
- Crystal symmetry
- Space group distributions
- Crystal system distributions
- Point group distributions
- Bravais lattices
- Hall symbols
- Symmetry operations
- Lattice parameters
- Density
- Cell volume
- Coordination environments
- Bond lengths
- Bond angles
- Radial distribution functions (RDF)
- Crystal fingerprints
- Structural similarity
- Novelty
- Dataset similarity
- Formation energy
- Energy above hull
- Band gap
- Stability metrics

The objective is to determine whether MatterGen generates physically meaningful, chemically realistic, and structurally diverse crystal structures compared with established experimental and computational materials databases.

---

## License

This project is licensed under the **MIT License**.
