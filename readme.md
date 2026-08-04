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
- Crystal system
- Structure validity

#### Statistics Output

- `basic_statistics.csv`
- `summary.csv`
- `failed_files.csv`
- Density histogram
- Volume histogram
- Number of atoms histogram
- Lattice parameter distributions
- Space group distribution
- Composition distribution

#### Statistics Example

```bash
python 01_basic_statistics.py \
    --input ../../data/MaterialsProject/Al-O \
    --output ../../results/MaterialsProject/Al-O \
    --workers -1
```

---

## Planned Analysis Pipeline

| Script | Analysis |
| ------ | -------- |
| 00_download_database.py | Download reference datasets |
| 01_basic_statistics.py | Basic structural statistics |
| 02_symmetry_analysis.py | Space groups and crystal systems |
| 03_coordination_analysis.py | CrystalNN and VoronoiNN coordination |
| 04_bond_analysis.py | Bond lengths, bond angles, nearest neighbors |
| 05_rdf_analysis.py | Radial distribution functions (RDF) |
| 06_chemical_analysis.py | Oxidation states and charge neutrality |
| 07_fingerprint_analysis.py | Crystal fingerprints |
| 08_novelty_analysis.py | Duplicate detection and structural similarity |
| 09_distribution_comparison.py | Statistical comparison between datasets |
| 10_stability_analysis.py | Formation energy, energy above hull, band gap and stability metrics |

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
│
├── MaterialsProject/
│   ├── Al-O/
│   ├── Fe-O/
│   ├── Ti-O/
│   └── Y-Al-O/
│
├── COD/
├── OQMD/
├── JARVIS/
└── Alexandria/

results/
├── MatterGen/
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
3. Run the same analysis pipeline on every dataset.
4. Compare statistical distributions.
5. Compare crystal symmetry distributions.
6. Compare local atomic environments.
7. Compare bond geometry and RDFs.
8. Evaluate novelty and structural diversity.
9. Compare stability-related properties where available.
10. Determine whether MatterGen generates realistic crystal structures.

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

The planned benchmark will compare:

- Chemical composition
- Crystal symmetry
- Space group distributions
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
