# MatterGen Benchmark Scripts

A collection of Python scripts for statistically benchmarking MatterGen-generated crystal structures against established crystal structure databases.

The benchmark is designed to compare structural, crystallographic, geometric, and statistical properties of generated materials in a reproducible and automated workflow. Every analysis produces publication-quality figures together with machine-readable CSV files for further statistical evaluation.

The project focuses exclusively on information that can be derived directly from crystal structures without requiring electronic structure calculations.

---

## Features

### Crystal Structure Analysis

- Analyze CIF, XYZ, and EXTXYZ crystal structures.
- Support thousands of crystal structures per dataset.
- Process multiple reference databases using a unified workflow.
- Automatically organize outputs into dedicated analysis directories.
- Support parallel processing on multicore systems.

### Statistical Analysis

- Analyze chemical composition.
- Analyze crystallographic symmetry.
- Analyze local coordination environments.
- Analyze bond geometry.
- Analyze radial distribution functions (RDF).
- Analyze crystal fingerprints.
- Compute statistical structural similarity.
- Perform principal component analysis (PCA).
- Export comprehensive CSV summaries.
- Generate publication-quality PDF and PNG figures.

### Benchmarking

The benchmark enables direct comparison between MatterGen-generated crystal structures and reference databases, including:

- Chemical composition distributions.
- Crystal systems.
- Space groups.
- Point groups.
- Hall symbols.
- Bravais lattices.
- Symmetry operations.
- Coordination environments.
- Bond-length distributions.
- Bond-angle distributions.
- Radial distribution functions.
- Crystal fingerprints.
- Structural similarity.
- Dataset diversity.

---

## Why Both RDF and Fingerprint Analysis?

Radial distribution functions (RDFs) and crystal fingerprints provide complementary information about crystal structures.

### Radial Distribution Functions

RDF analysis describes the statistical distribution of interatomic distances within a dataset. It answers questions such as:

- How far apart are neighboring atoms?
- How similar are average atomic distances?
- Are generated structures statistically consistent with reference datasets?

RDFs characterize average geometric behavior but do not directly encode local coordination geometry.

### Crystal Fingerprints

Crystal fingerprints encode local atomic environments into fixed-length numerical descriptors. These descriptors capture information about:

- Local coordination.
- Bond topology.
- Neighbor chemistry.
- Local geometric environments.

Fingerprints enable quantitative comparison of complete crystal structures using similarity metrics such as cosine similarity.

Together, RDF and fingerprint analysis provide complementary views of structural similarity.

---

## Table of Contents

- [MatterGen Benchmark Scripts](#mattergen-benchmark-scripts)
  - [Features](#features)
    - [Crystal Structure Analysis](#crystal-structure-analysis)
    - [Statistical Analysis](#statistical-analysis)
    - [Benchmarking](#benchmarking)
  - [Why Both RDF and Fingerprint Analysis?](#why-both-rdf-and-fingerprint-analysis)
    - [Radial Distribution Functions](#radial-distribution-functions)
    - [Crystal Fingerprints](#crystal-fingerprints)
  - [Table of Contents](#table-of-contents)
  - [Requirements](#requirements)
    - [Python](#python)
    - [Required Packages](#required-packages)
  - [Installation](#installation)
  - [Project Structure](#project-structure)
  - [Example Dataset Directory Structure](#example-dataset-directory-structure)
  - [Output Directory Structure](#output-directory-structure)
  - [Supported Databases](#supported-databases)
  - [Supported File Formats](#supported-file-formats)
  - [Command-Line Interface](#command-line-interface)
    - [Arguments](#arguments)
  - [General Workflow](#general-workflow)
  - [Analysis Pipeline](#analysis-pipeline)
- [Structural Benchmark Methodology](#structural-benchmark-methodology)
  - [1. Chemical Composition](#1-chemical-composition)
  - [2. Crystallographic Symmetry](#2-crystallographic-symmetry)
  - [3. Local Atomic Geometry](#3-local-atomic-geometry)
  - [4. Statistical Structure Descriptors](#4-statistical-structure-descriptors)
    - [Radial Distribution Functions (RDF)](#radial-distribution-functions-rdf)
    - [Crystal Fingerprints](#crystal-fingerprints-1)
      - [CrystalNN Fingerprint](#crystalnn-fingerprint)
      - [OPSite Fingerprint](#opsite-fingerprint)
      - [Voronoi Fingerprint](#voronoi-fingerprint)
  - [Structural Similarity](#structural-similarity)
  - [Principal Component Analysis (PCA)](#principal-component-analysis-pca)
  - [Computational Complexity](#computational-complexity)
  - [Parallel Processing](#parallel-processing)
- [01\_composition\_analysis.py](#01_composition_analysispy)
  - [Purpose](#purpose)
  - [Computed Properties](#computed-properties)
  - [Generated Files](#generated-files)
  - [Example](#example)
  - [Interpretation](#interpretation)
- [02\_symmetry\_analysis.py](#02_symmetry_analysispy)
  - [Purpose](#purpose-1)
  - [Computed Properties](#computed-properties-1)
  - [Generated Files](#generated-files-1)
  - [Example](#example-1)
  - [Interpretation](#interpretation-1)
  - [Notes](#notes)
- [03\_coordination\_analysis.py](#03_coordination_analysispy)
  - [Purpose](#purpose-2)
  - [Neighbor-Finding Algorithms](#neighbor-finding-algorithms)
    - [CrystalNN](#crystalnn)
    - [VoronoiNN](#voronoinn)
  - [Computed Properties](#computed-properties-2)
  - [Generated Files](#generated-files-2)
  - [Example](#example-2)
  - [Interpretation](#interpretation-2)
- [04\_bond\_analysis.py](#04_bond_analysispy)
  - [Purpose](#purpose-3)
  - [Computed Properties](#computed-properties-3)
  - [Generated Files](#generated-files-3)
  - [Example](#example-3)
  - [Interpretation](#interpretation-3)
  - [Relationship Between Local Analyses](#relationship-between-local-analyses)
- [05\_rdf\_analysis.py](#05_rdf_analysispy)
  - [Purpose](#purpose-4)
  - [Computed Properties](#computed-properties-4)
  - [Generated Files](#generated-files-4)
  - [Example](#example-4)
  - [Interpretation](#interpretation-4)
- [06\_fingerprint\_analysis.py](#06_fingerprint_analysispy)
  - [Purpose](#purpose-5)
  - [Fingerprint Types](#fingerprint-types)
    - [CrystalNN Fingerprint](#crystalnn-fingerprint-1)
    - [OPSite Fingerprint](#opsite-fingerprint-1)
    - [Voronoi Fingerprint](#voronoi-fingerprint-1)
  - [Structural Similarity](#structural-similarity-1)
    - [Cosine Similarity](#cosine-similarity)
  - [Principal Component Analysis](#principal-component-analysis)
  - [Computed Properties](#computed-properties-5)
  - [Generated Files](#generated-files-5)
  - [Example](#example-5)
  - [Interpretation](#interpretation-5)
  - [RDF vs. Fingerprint Analysis](#rdf-vs-fingerprint-analysis)
- [Output Structure](#output-structure)
- [Typical Benchmark Workflow](#typical-benchmark-workflow)
  - [Step 1](#step-1)
  - [Step 2](#step-2)
  - [Step 3](#step-3)
  - [Step 4](#step-4)
  - [Step 5](#step-5)
  - [Step 6](#step-6)
- [Interpreting the Results](#interpreting-the-results)
- [Performance](#performance)
- [Limitations](#limitations)
- [Troubleshooting](#troubleshooting)
  - [CrystalNN Warning](#crystalnn-warning)
  - [XYZ Files](#xyz-files)
  - [Invalid Structures](#invalid-structures)
  - [Parallel Processing](#parallel-processing-1)
- [References](#references)
  - [Core Libraries](#core-libraries)
  - [Materials Science](#materials-science)
  - [Machine Learning](#machine-learning)
- [License](#license)

---

## Requirements

### Python

- Python 3.11 or newer

### Required Packages

| Package | Purpose |
| ---------- | --------- |
| numpy | Numerical computations |
| scipy | Scientific computing |
| pandas | Data analysis |
| matplotlib | Plot generation |
| seaborn | Statistical visualization |
| pymatgen | Crystal structure analysis |
| pymatviz | Radial distribution functions |
| matminer | Crystal fingerprints |
| ase | Reading crystal structure files |
| spglib | Symmetry analysis |
| scikit-learn | PCA and similarity analysis |
| joblib | Parallel processing |
| tqdm | Progress bars |

Install all dependencies using:

```bash
pip install numpy scipy pandas matplotlib seaborn pymatgen pymatviz matminer ase spglib scikit-learn joblib tqdm
```

---

## Installation

Clone the repository:

```bash
git clone <repository-url>
cd MatterGen-Scripts
```

Install all required dependencies:

```bash
pip install -r requirements.txt
```

If a requirements file is unavailable, install the packages listed above manually.

## Project Structure

```text
MatterGen-Scripts/
├── 01_composition_analysis.py
├── 02_symmetry_analysis.py
├── 03_coordination_analysis.py
├── 04_bond_analysis.py
├── 05_rdf_analysis.py
├── 06_fingerprint_analysis.py
├── requirements.txt
├── README.md
├── data/
├── results/
└── .vscode/
```

Each script performs one independent statistical analysis and stores its results in a dedicated subdirectory of the selected output directory.

The scripts are intentionally independent, allowing individual analyses to be executed without running the complete benchmark.

---

## Example Dataset Directory Structure

The benchmark expects the following directory layout.

```text
data/
├── MatterGen/
│   ├── Al-O/
│   │   └── generated_crystals_cif/
│   ├── Fe-O/
│   │   ├── generated_crystals_cif/
│   ├── Ti-O/
│   │   └── generated_crystals_cif/
│   └── Y-Al-O/
│       └── generated_crystals_cif/
│
├── MaterialsProject/
│   ├── Al-O/
│   ├── Fe-O/
│   ├── Ti-O/
│   └── Y-Al-O/
│
├── OQMD/
│   ├── Al-O/
│   ├── Fe-O/
│   ├── Ti-O/
│   └── Y-Al-O/
│
├── JARVIS/
│   ├── Al-O/
│   ├── Fe-O/
│   ├── Ti-O/
│   └── Y-Al-O/
│
├── Alexandria/
│   ├── Al-O/
│   ├── Fe-O/
│   ├── Ti-O/
│   └── Y-Al-O/
│
└── COD/
    ├── Al-O/
    ├── Fe-O/
    ├── Ti-O/
    └── Y-Al-O/
```

Each chemical system is analyzed independently.

---

## Output Directory Structure

Each analysis script creates its own output directory.

```text
results/
└── <database>/
    └── <chemical-system>/
        ├── composition/
        ├── symmetry/
        ├── coordination/
        ├── bond/
        ├── rdf/
        └── fingerprint/
```

Each analysis directory contains:

- CSV files
- PDF figures
- PNG figures
- summary tables
- log files for failed structures (when applicable)

The analyses are completely independent and never overwrite the outputs of other scripts.

---

## Supported Databases

The benchmark currently supports the following datasets.

| Database | Purpose |
| ---------- | --------- |
| MatterGen | Generated crystal structures |
| Materials Project | DFT-computed reference structures |
| OQMD | Large computational materials database |
| JARVIS | DFT materials database |
| Alexandria | Computational materials database |
| COD | Experimental crystal structures |

Additional databases can be benchmarked by providing structures in one of the supported file formats.

---

## Supported File Formats

The scripts automatically detect the following structure formats.

| Extension | Description |
| ----------- | ------------- |
| `.cif` | Crystallographic Information File |
| `.xyz` | XYZ atomic coordinates |
| `.extxyz` | Extended XYZ with lattice information |

The benchmark recursively searches input directories when the `--recursive` option is specified.

---

## Command-Line Interface

Every analysis script follows the same command-line interface.

```text
python <script>.py \
    --input <input_directory> \
    --output <output_directory> \
    --workers -1 \
    --recursive
```

### Arguments

| Argument | Description |
| ---------- | ------------- |
| `--input` | Directory containing crystal structures |
| `--output` | Output directory |
| `--workers` | Number of CPU cores (`-1` uses all available cores) |
| `--recursive` | Search subdirectories recursively |

This consistent interface allows all benchmark scripts to be executed interchangeably using the same workflow.

---

## General Workflow

Every analysis follows the same sequence of operations.

```text
Input structures
        │
        ▼
Read crystal structures
        │
        ▼
Validate structures
        │
        ▼
Compute structural descriptors
        │
        ▼
Statistical analysis
        │
        ▼
Generate CSV files
        │
        ▼
Generate publication-quality figures
        │
        ▼
Write summary statistics
```

This workflow ensures consistent benchmarking across all implemented analyses.

## Analysis Pipeline

The benchmark is organized into independent analysis modules. Each module focuses on a specific structural aspect and produces dedicated CSV files and publication-quality figures.

| Script | Analysis | Output Directory |
| ---------- | ---------- | ------------------ |
| `01_composition_analysis.py` | Chemical composition statistics | `<chemical-system>/` |
| `02_symmetry_analysis.py` | Crystallographic symmetry | `symmetry/` |
| `03_coordination_analysis.py` | Local coordination environments | `coordination/` |
| `04_bond_analysis.py` | Bond lengths and bond angles | `bond/` |
| `05_rdf_analysis.py` | Radial distribution functions | `rdf/` |
| `06_fingerprint_analysis.py` | Crystal fingerprints and structural similarity | `fingerprint/` |

All scripts share the same command-line interface and produce consistent output formats, allowing them to be combined into a complete benchmarking workflow.

---

# Structural Benchmark Methodology

The benchmark compares MatterGen-generated crystal structures with reference materials databases using complementary structural descriptors.

The implemented analyses can be grouped into four categories.

## 1. Chemical Composition

Composition analysis evaluates the elemental makeup of generated structures.

Typical quantities include:

- Chemical formulas
- Reduced formulas
- Element frequencies
- Stoichiometric distributions
- Number of atoms per structure

Composition analysis answers questions such as:

- Are the generated compositions chemically reasonable?
- Are certain elements overrepresented?
- Does the composition distribution match the reference dataset?

---

## 2. Crystallographic Symmetry

Symmetry analysis evaluates the crystallographic properties of crystal structures.

Computed properties include:

- Crystal system
- Space group
- Point group
- Hall symbol
- Bravais lattice
- Number of symmetry operations

Symmetry analysis answers questions such as:

- Which crystal systems dominate the dataset?
- Does MatterGen reproduce realistic symmetry distributions?
- Are low-symmetry structures overgenerated?

---

## 3. Local Atomic Geometry

Local geometry characterizes the arrangement of neighboring atoms.

This category includes:

- Coordination analysis
- Bond-length analysis
- Bond-angle analysis

These analyses describe the local atomic environment independent of the overall crystal symmetry.

---

## 4. Statistical Structure Descriptors

The benchmark implements two complementary statistical descriptors.

### Radial Distribution Functions (RDF)

Radial distribution functions measure the probability of finding atoms at a given distance from another atom.

RDFs describe:

- Average bond distances
- Long-range order
- Dataset-wide structural statistics

RDFs are particularly useful for comparing the overall geometric distributions of different datasets.

---

### Crystal Fingerprints

Crystal fingerprints encode local atomic environments into numerical feature vectors.

The benchmark computes three fingerprint types.

#### CrystalNN Fingerprint

Based on CrystalNN nearest-neighbor analysis.

Captures:

- Coordination geometry
- Local bonding environment
- Neighbor chemistry

---

#### OPSite Fingerprint

Based on order parameters describing ideal coordination geometries.

Captures similarity to environments such as:

- Tetrahedral
- Octahedral
- Trigonal planar
- Square planar
- Cubic
- Icosahedral

---

#### Voronoi Fingerprint

Based on Voronoi tessellation.

Captures:

- Voronoi coordination
- Neighbor topology
- Local geometric environments

Unlike CrystalNN fingerprints, Voronoi fingerprints do not require predefined bonding criteria.

---

## Structural Similarity

Crystal fingerprints enable direct comparison between complete crystal structures.

The benchmark computes cosine similarity between fingerprint vectors.

Interpretation:

| Cosine Similarity | Interpretation |
| ------------------ | --------------- |
| 1.00 | Nearly identical local environments |
| 0.90–0.99 | Very similar structures |
| 0.70–0.90 | Similar local coordination |
| 0.50–0.70 | Moderately similar |
| < 0.50 | Structurally different |

Higher similarity indicates more comparable local atomic environments.

---

## Principal Component Analysis (PCA)

Fingerprint vectors typically contain hundreds of numerical features.

Principal Component Analysis (PCA) projects these high-dimensional fingerprints into two dimensions for visualization.

PCA plots reveal:

- Clusters of similar structures
- Dataset diversity
- Outliers
- Structural trends

Nearby points in the PCA plot generally correspond to structurally similar crystals.

---

## Computational Complexity

Approximate computational complexity of each analysis.

| Script | Complexity |
| --------- | ----------- |
| Composition | O(N) |
| Symmetry | O(N) |
| Coordination | O(N × neighbors) |
| Bond analysis | O(N × bonds) |
| RDF | O(N × atoms²) |
| Fingerprints | O(N × atoms × neighbors) |

where **N** is the number of crystal structures.

The fingerprint analysis is generally the most computationally intensive due to the calculation of local atomic environments for every atomic site.

---

## Parallel Processing

All benchmark scripts support parallel execution using the `--workers` option.

Example:

```bash
python 06_fingerprint_analysis.py \
    --input data/MatterGen/Al-O/generated_crystals_cif \
    --output results/MatterGen/Al-O \
    --workers -1
```

Using

```text
--workers -1
```

automatically utilizes all available CPU cores.

For large datasets containing several thousand crystal structures, parallel execution can significantly reduce processing time.

# 01_composition_analysis.py

Analyzes the chemical composition of crystal structure datasets and generates descriptive statistics for elemental distributions and stoichiometry.

## Purpose

Composition analysis verifies that generated crystal structures exhibit realistic chemical compositions and enables comparison with experimental and computational reference databases.

The analysis focuses exclusively on information contained within the crystal structures and does not require electronic structure calculations.

---

## Computed Properties

For each crystal structure, the script determines:

- Chemical formula
- Reduced chemical formula
- Number of unique elements
- Number of atoms
- Elemental composition
- Stoichiometric ratios

Dataset-wide statistics include:

- Element frequency
- Formula frequency
- Reduced formula frequency
- Distribution of atoms per structure
- Distribution of elements per structure

---

## Generated Files

```text
composition/
├── composition_analysis.csv
├── formula_distribution.csv
├── reduced_formula_distribution.csv
├── element_distribution.csv
├── summary.csv
├── failed_files.csv
├── formula_distribution.pdf
├── formula_distribution.png
├── element_distribution.pdf
├── element_distribution.png
├── atoms_per_structure.pdf
└── atoms_per_structure.png
```

---

## Example

```bash
python 01_composition_analysis.py \
    --input ../../data/MaterialsProject/Al-O \
    --output ../../results/MaterialsProject/Al-O \
    --workers -1
```

---

## Interpretation

Useful questions answered by this analysis include:

- Are the generated compounds chemically reasonable?
- Does the elemental distribution resemble the reference dataset?
- Are some compositions overrepresented?
- Does the generator produce chemically diverse structures?

---

# 02_symmetry_analysis.py

Analyzes crystallographic symmetry using `spglib` and `pymatgen`.

## Purpose

Symmetry analysis compares the crystallographic characteristics of generated crystal structures with established materials databases.

Crystal symmetry is one of the most important indicators of structural realism and diversity.

---

## Computed Properties

For every crystal structure:

- Crystal system
- Space group
- Space group number
- Point group
- Hall symbol
- Bravais lattice
- Number of symmetry operations

Dataset-wide statistics include:

- Crystal-system distribution
- Space-group distribution
- Point-group distribution
- Hall-symbol distribution
- Bravais-lattice distribution

---

## Generated Files

```text
symmetry/
├── symmetry_analysis.csv
├── crystal_system_distribution.csv
├── space_group_distribution.csv
├── point_group_distribution.csv
├── hall_symbol_distribution.csv
├── bravais_lattice_distribution.csv
├── summary.csv
├── failed_files.csv
├── crystal_system_distribution.pdf
├── crystal_system_distribution.png
├── space_group_distribution.pdf
├── space_group_distribution.png
├── point_group_distribution.pdf
├── point_group_distribution.png
├── bravais_lattice_distribution.pdf
└── bravais_lattice_distribution.png
```

---

## Example

```bash
python 02_symmetry_analysis.py \
    --input ../../data/MaterialsProject/Al-O \
    --output ../../results/MaterialsProject/Al-O \
    --workers -1
```

---

## Interpretation

This analysis answers questions such as:

- Does MatterGen reproduce realistic crystal-system distributions?
- Which space groups are most common?
- Are generated structures predominantly high- or low-symmetry?
- How does the symmetry distribution compare with reference databases?

Large deviations from reference distributions may indicate biases in the structure generation process.

---

## Notes

Crystal symmetry is determined directly from the atomic positions and lattice parameters contained in each structure.

The reported symmetry may depend on the numerical tolerance used by the underlying symmetry detection algorithms.

# 03_coordination_analysis.py

Analyzes local atomic coordination environments using multiple nearest-neighbor algorithms.

## Purpose

Coordination analysis characterizes the local atomic environments of crystal structures. It determines how many neighboring atoms surround each atom and compares the statistical coordination distributions between MatterGen-generated structures and reference materials databases.

Unlike symmetry analysis, coordination analysis focuses on local atomic environments rather than the global crystal structure.

---

## Neighbor-Finding Algorithms

The benchmark computes coordination numbers using two independent methods.

### CrystalNN

CrystalNN determines neighboring atoms using chemically informed bonding heuristics implemented in `pymatgen`.

Advantages:

- Chemically meaningful neighbors
- Robust for most inorganic crystals
- Considers atomic radii and bonding environments

---

### VoronoiNN

VoronoiNN determines neighbors using Voronoi tessellation.

Advantages:

- Purely geometry based
- Independent of bonding assumptions
- Useful for comparing diverse crystal structures

---

## Computed Properties

For every crystal structure:

- Average CrystalNN coordination number
- Average VoronoiNN coordination number
- Minimum coordination number
- Maximum coordination number
- Coordination-number distribution

Dataset-wide statistics include:

- Mean coordination number
- Median coordination number
- Standard deviation
- Coordination histograms
- Coordination frequency distributions

---

## Generated Files

```text
coordination/
├── coordination_analysis.csv
├── crystalnn_distribution.csv
├── voronoi_distribution.csv
├── summary.csv
├── failed_files.csv
├── crystalnn_histogram.pdf
├── crystalnn_histogram.png
├── voronoi_histogram.pdf
├── voronoi_histogram.png
├── coordination_boxplot.pdf
└── coordination_boxplot.png
```

---

## Example

```bash
python 03_coordination_analysis.py \
    --input ../../data/MaterialsProject/Al-O \
    --output ../../results/MaterialsProject/Al-O \
    --workers -1
```

---

## Interpretation

Coordination analysis answers questions such as:

- Are generated coordination environments chemically realistic?
- Does MatterGen reproduce the coordination statistics of known materials?
- Which coordination numbers dominate the dataset?
- Are under- or over-coordinated environments common?

Large deviations from reference datasets may indicate unrealistic local atomic environments.

---

# 04_bond_analysis.py

Analyzes bond lengths and bond angles of crystal structures.

## Purpose

Bond geometry provides a direct measure of local structural realism.

This analysis evaluates statistical distributions of bond lengths and bond angles across complete datasets.

Bond geometry complements coordination analysis by describing not only the number of neighbors but also their geometric arrangement.

---

## Computed Properties

For every crystal structure:

- Bond lengths
- Bond-angle values
- Mean bond length
- Mean bond angle
- Minimum bond length
- Maximum bond length
- Minimum bond angle
- Maximum bond angle

Dataset-wide statistics include:

- Bond-length distribution
- Bond-angle distribution
- Mean bond-length statistics
- Mean bond-angle statistics

---

## Generated Files

```text
bond/
├── bond_analysis.csv
├── bond_lengths.csv
├── bond_angles.csv
├── summary.csv
├── failed_files.csv
├── bond_length_histogram.pdf
├── bond_length_histogram.png
├── bond_angle_histogram.pdf
├── bond_angle_histogram.png
├── bond_length_boxplot.pdf
├── bond_length_boxplot.png
├── bond_angle_boxplot.pdf
└── bond_angle_boxplot.png
```

---

## Example

```bash
python 04_bond_analysis.py \
    --input ../../data/MaterialsProject/Al-O \
    --output ../../results/MaterialsProject/Al-O \
    --workers -1
```

---

## Interpretation

Bond analysis helps answer questions such as:

- Are bond lengths physically reasonable?
- Does MatterGen reproduce realistic bond-angle distributions?
- Are generated structures excessively distorted?
- Does the local geometry resemble known crystal structures?

Well-aligned bond-length and bond-angle distributions indicate that the generated structures closely resemble experimentally observed or computationally predicted materials.

---

## Relationship Between Local Analyses

The three local-structure analyses complement each other.

| Analysis | Measures |
| ---------- | ---------- |
| Coordination | Number of neighboring atoms |
| Bond analysis | Local geometry and bond distances |
| Fingerprint analysis | Numerical description of complete local environments |

Together, these analyses provide a comprehensive characterization of local atomic structure while remaining independent of global crystallographic symmetry.

# 05_rdf_analysis.py

Analyzes radial distribution functions (RDFs) to compare the statistical atomic-distance distributions of crystal structure datasets.

## Purpose

Radial distribution functions describe how atomic density varies as a function of distance from a reference atom.

Unlike bond analysis, which considers only bonded neighbors, RDF analysis captures short-, medium-, and long-range structural order.

RDFs provide a statistical fingerprint of an entire dataset and are widely used to compare generated materials with experimental and computational reference databases.

---

## Computed Properties

For every crystal structure:

- Global radial distribution function
- RDF distance bins
- RDF intensity values

Dataset-wide statistics include:

- Mean RDF
- Median RDF
- Standard deviation
- Minimum RDF
- Maximum RDF

---

## Generated Files

```text
rdf/
├── rdf_analysis.csv
├── rdf_values.csv
├── mean_rdf.csv
├── summary.csv
├── failed_files.csv
├── rdf_plot.pdf
├── rdf_plot.png
├── rdf_heatmap.pdf
├── rdf_heatmap.png
├── rdf_boxplot.pdf
└── rdf_boxplot.png
```

---

## Example

```bash
python 05_rdf_analysis.py \
    --input ../../data/MaterialsProject/Al-O \
    --output ../../results/MaterialsProject/Al-O \
    --workers -1
```

---

## Interpretation

RDF analysis answers questions such as:

- Are average atomic distances realistic?
- Does MatterGen reproduce the overall structural statistics of known materials?
- Are long-range ordering characteristics preserved?
- Which datasets exhibit similar geometric distributions?

Because RDFs average over all atomic pairs, they describe dataset-level structural behavior rather than individual coordination environments.

---

# 06_fingerprint_analysis.py

Analyzes crystal fingerprints and structural similarity using local-environment descriptors from **matminer**.

## Purpose

Crystal fingerprints encode local atomic environments as numerical feature vectors.

Unlike RDF analysis, which measures statistical atomic-distance distributions, fingerprint analysis captures local coordination geometry and enables direct quantitative comparison between complete crystal structures.

The benchmark computes three complementary fingerprint types to characterize structural diversity and similarity.

---

## Fingerprint Types

### CrystalNN Fingerprint

Based on chemically informed nearest-neighbor analysis.

Captures:

- Local coordination
- Neighbor chemistry
- Bond topology
- Coordination geometry

---

### OPSite Fingerprint

Based on order parameters describing ideal coordination environments.

Measures similarity to environments such as:

- Tetrahedral
- Octahedral
- Trigonal planar
- Square planar
- Cubic
- Body-centered cubic
- Icosahedral

---

### Voronoi Fingerprint

Based on Voronoi tessellation of the crystal structure.

Captures:

- Geometric neighbor topology
- Voronoi coordination
- Local packing environment

Unlike CrystalNN, Voronoi fingerprints do not depend on bonding heuristics.

---

## Structural Similarity

The benchmark compares fingerprint vectors using cosine similarity.

### Cosine Similarity

Cosine similarity measures the angular similarity between two fingerprint vectors.

| Similarity | Interpretation |
| -----------:| --------------- |
| 1.00 | Nearly identical structures |
| 0.90–0.99 | Very similar local environments |
| 0.70–0.90 | Similar coordination environments |
| 0.50–0.70 | Moderately similar |
| < 0.50 | Structurally different |

Higher cosine similarity indicates more similar local atomic environments.

---

## Principal Component Analysis

Principal Component Analysis (PCA) projects high-dimensional fingerprint vectors into two dimensions.

The resulting PCA plots help identify:

- Structural clusters
- Dataset diversity
- Outliers
- Relationships between structures

Nearby points generally correspond to structurally similar crystals.

---

## Computed Properties

For every crystal structure:

- CrystalNN fingerprint
- OPSite fingerprint
- Voronoi fingerprint
- CrystalNN fingerprint norm
- OPSite fingerprint norm
- Voronoi fingerprint norm
- Mean cosine similarity
- Nearest-neighbor similarity
- Principal component coordinates

Dataset-wide statistics include:

- Mean fingerprint norms
- Median fingerprint norms
- Fingerprint similarity distribution
- Pairwise similarity matrix
- PCA embedding
- Dataset diversity statistics

---

## Generated Files

```text
fingerprint/
├── fingerprint_analysis.csv
├── crystalnn_fingerprints.csv
├── opsite_fingerprints.csv
├── voronoi_fingerprints.csv
├── similarity_matrix.csv
├── pca_coordinates.csv
├── summary.csv
├── failed_files.csv
├── fingerprint_similarity_histogram.pdf
├── fingerprint_similarity_histogram.png
├── fingerprint_heatmap.pdf
├── fingerprint_heatmap.png
├── fingerprint_pca.pdf
├── fingerprint_pca.png
├── fingerprint_boxplot.pdf
└── fingerprint_boxplot.png
```

---

## Example

```bash
python 06_fingerprint_analysis.py \
    --input ../../data/MaterialsProject/Al-O \
    --output ../../results/MaterialsProject/Al-O \
    --workers -1
```

---

## Interpretation

Fingerprint analysis answers questions such as:

- Are MatterGen structures locally similar to reference materials?
- How diverse are the generated structures?
- Do clusters of similar crystal structures exist?
- Are duplicate or nearly identical structures generated?
- Which structures are statistical outliers?

Unlike RDF analysis, fingerprint analysis compares complete local atomic environments rather than average interatomic distance distributions.

---

## RDF vs. Fingerprint Analysis

Although both analyses describe structural similarity, they capture different information.

| RDF Analysis | Fingerprint Analysis |
| -------------- | ---------------------- |
| Average atomic-distance distribution | Local atomic environments |
| Long-range statistical order | Local coordination geometry |
| Dataset-level descriptor | Structure-level descriptor |
| Continuous distance distribution | High-dimensional feature vector |
| Compares average geometry | Compares complete local environments |

Using both analyses provides a more comprehensive assessment of structural realism than either method alone.

# Output Structure

Each analysis creates its own dedicated output directory.

```text
results/
└── <database>/
    └── <chemical-system>/
        ├── composition/
        ├── symmetry/
        ├── coordination/
        ├── bond/
        ├── rdf/
        └── fingerprint/
```

Each analysis directory contains one or more of the following:

- Analysis CSV files
- Statistical summary tables
- Publication-quality PDF figures
- High-resolution PNG figures
- Failed file reports

The generated CSV files are intended for further statistical analysis, while the figures provide immediate visual comparison between datasets.

---

# Typical Benchmark Workflow

A complete benchmark typically consists of the following steps.

## Step 1

Analyze chemical composition.

```text
01_composition_analysis.py
```

Verify that generated structures contain realistic compositions.

---

## Step 2

Analyze crystallographic symmetry.

```text
02_symmetry_analysis.py
```

Compare crystal systems, space groups, and symmetry statistics.

---

## Step 3

Analyze local coordination.

```text
03_coordination_analysis.py
```

Compare nearest-neighbor environments.

---

## Step 4

Analyze bond geometry.

```text
04_bond_analysis.py
```

Compare bond-length and bond-angle distributions.

---

## Step 5

Analyze radial distribution functions.

```text
05_rdf_analysis.py
```

Compare average structural statistics.

---

## Step 6

Analyze crystal fingerprints.

```text
06_fingerprint_analysis.py
```

Compare local atomic environments and structural similarity.

---

# Interpreting the Results

Each analysis measures a different structural property.

| Analysis | Interpretation |
| ---------- | ---------------- |
| Composition | Chemical realism |
| Symmetry | Crystallographic realism |
| Coordination | Local neighbor environments |
| Bond analysis | Local geometry |
| RDF | Statistical distance distributions |
| Fingerprints | Structural similarity |

No individual analysis is sufficient on its own.

The benchmark is intended to combine all analyses to obtain a comprehensive evaluation of generated crystal structures.

---

# Performance

Approximate relative runtime.

| Analysis | Runtime |
| ---------- | --------- |
| Composition | Very Fast |
| Symmetry | Fast |
| Coordination | Moderate |
| Bond Analysis | Moderate |
| RDF | Slow |
| Fingerprint | Slow |

Fingerprint analysis is generally the most computationally expensive because local descriptors are computed for every atomic site.

Using

```text
--workers -1
```

is recommended for large datasets.

---

# Limitations

The benchmark analyzes structural information that is directly available from crystal structure files.

The following quantities **cannot** be computed reliably from CIF, XYZ, or EXTXYZ files alone:

- Formation energy
- Energy above hull
- Electronic band gap
- Elastic constants
- Magnetic properties
- Phonon properties
- Electronic density of states

These quantities require additional electronic structure calculations (for example, density functional theory) or values obtained from reference materials databases.

---

# Troubleshooting

## CrystalNN Warning

Example:

```text
CrystalNN: cannot locate an appropriate radius...
```

This warning indicates that ionic radii are unavailable for one or more atoms.

`CrystalNN` automatically falls back to covalent or atomic radii.

The warning does **not** invalidate the fingerprint analysis.

---

## XYZ Files

Standard XYZ files usually do not contain lattice vectors.

Fingerprint, symmetry, and RDF analyses require periodic crystal structures.

Whenever possible, use:

- CIF
- EXTXYZ

instead of plain XYZ files.

---

## Invalid Structures

Malformed crystal structures are skipped automatically.

Details are written to:

```text
failed_files.csv
```

---

## Parallel Processing

If parallel execution causes excessive memory usage, reduce the number of workers.

Example:

```bash
--workers 8
```

instead of

```bash
--workers -1
```

---

# References

If you use this benchmark in scientific work, please cite the software packages used by the analysis scripts.

## Core Libraries

- NumPy
- SciPy
- pandas
- matplotlib
- seaborn

## Materials Science

- pymatgen
- pymatviz
- matminer
- spglib
- ASE

## Machine Learning

- scikit-learn

Please also cite the corresponding materials databases used for benchmarking, such as:

- MatterGen
- Materials Project
- OQMD
- JARVIS
- Alexandria
- Crystallography Open Database (COD)

---

# License

This project is intended for academic research and benchmarking of crystal structure generation models.

Please refer to the repository license for usage and redistribution terms.