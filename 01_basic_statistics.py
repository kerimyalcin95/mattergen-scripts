#!/usr/bin/env python3
"""
01_basic_statistics.py

Compute basic statistics for a dataset of crystal structures.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes:
    * Chemical formula
    * Number of atoms
    * Lattice parameters
    * Cell volume
    * Density
    * Space group
- Saves:
    * CSV with one row per structure
    * Summary CSV
    * Publication-quality PDF figures
- Prints summary to terminal

Supported libraries:
    pymatgen
    ase
    numpy
    pandas
    scipy
    matplotlib
    seaborn
    joblib
    tqdm
"""

from __future__ import annotations

import argparse
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

from tqdm import tqdm
from joblib import Parallel, delayed

import matplotlib.pyplot as plt
import seaborn as sns

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ase.io import read as ase_read

# -----------------------------------------------------------------------------

SUPPORTED_EXTENSIONS = {
    ".cif",
    ".xyz",
    ".extxyz",
}

# -----------------------------------------------------------------------------

sns.set_theme(style="whitegrid")

plt.rcParams["figure.figsize"] = (8, 5)
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["axes.titlesize"] = 14
plt.rcParams["axes.labelsize"] = 12

# -----------------------------------------------------------------------------


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description="Compute basic statistics of crystal datasets."
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input folder containing crystal structures.",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output directory.",
    )

    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search subdirectories recursively.",
    )

    parser.add_argument(
        "--workers",
        default=-1,
        type=int,
        help="-1 uses all CPU cores.",
    )

    return parser.parse_args()


# -----------------------------------------------------------------------------


def find_structure_files(
    folder: Path,
    recursive: bool,
):
    """
    Find supported structure files.

    Parameters
    ----------
    folder
        Input directory.

    recursive
        Search recursively.

    Returns
    -------
    list[Path]
    """

    if recursive:
        iterator = folder.rglob("*")
    else:
        iterator = folder.glob("*")

    files = []

    for path in iterator:
        if path.suffix.lower() in SUPPORTED_EXTENSIONS:
            files.append(path)

    files.sort()

    return files


# -----------------------------------------------------------------------------


def load_structure(path: Path):
    """
    Load a crystal structure using pymatgen or ASE.

    Parameters
    ----------
    path

    Returns
    -------
    pymatgen Structure
    """

    suffix = path.suffix.lower()

    if suffix == ".cif":
        return Structure.from_file(path)

    atoms = ase_read(path)

    return AseAtomsAdaptor.get_structure(atoms)


# -----------------------------------------------------------------------------


def ensure_output_directory(
    output_dir: Path,
):
    """
    Create output directory if necessary.
    """

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


# -----------------------------------------------------------------------------


def safe_spacegroup(
    structure: Structure,
):
    """
    Determine the space group.

    Returns
    -------
    str
    """

    try:
        analyzer = SpacegroupAnalyzer(
            structure,
            symprec=0.1,
        )

        return analyzer.get_space_group_symbol()

    except Exception:
        return "Unknown"


# -----------------------------------------------------------------------------


def safe_density(
    structure: Structure,
):
    """
    Compute density.

    Returns
    -------
    float
    """

    try:
        return float(structure.density)

    except Exception:
        return np.nan


# -----------------------------------------------------------------------------


def safe_volume(
    structure: Structure,
):
    """
    Compute cell volume.
    """

    try:
        return float(structure.volume)

    except Exception:
        return np.nan

# -----------------------------------------------------------------------------


def analyze_structure(path: Path):
    """
    Analyze a single structure.

    Parameters
    ----------
    path : Path

    Returns
    -------
    dict
        Dictionary containing all computed properties.
    """

    try:
        structure = load_structure(path)

        lattice = structure.lattice

        result = {
            "file": path.name,
            "formula": structure.composition.formula,
            "reduced_formula": structure.composition.reduced_formula,
            "num_atoms": len(structure),
            "num_elements": len(structure.composition.elements),
            "a": lattice.a,
            "b": lattice.b,
            "c": lattice.c,
            "alpha": lattice.alpha,
            "beta": lattice.beta,
            "gamma": lattice.gamma,
            "volume": safe_volume(structure),
            "density": safe_density(structure),
            "space_group": safe_spacegroup(structure),
            "valid": True,
        }

        return result

    except Exception as exc:

        return {
            "file": path.name,
            "formula": None,
            "reduced_formula": None,
            "num_atoms": np.nan,
            "num_elements": np.nan,
            "a": np.nan,
            "b": np.nan,
            "c": np.nan,
            "alpha": np.nan,
            "beta": np.nan,
            "gamma": np.nan,
            "volume": np.nan,
            "density": np.nan,
            "space_group": "Invalid",
            "valid": False,
            "error": str(exc),
        }


# -----------------------------------------------------------------------------


def analyze_dataset(
    files,
    workers,
):
    """
    Analyze an entire dataset in parallel.

    Parameters
    ----------
    files
        List of structure files.

    workers
        Number of worker processes.

    Returns
    -------
    pandas.DataFrame
    """

    iterator = tqdm(
        files,
        desc="Analyzing structures",
    )

    rows = Parallel(
        n_jobs=workers,
    )(
        delayed(analyze_structure)(f)
        for f in iterator
    )

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------


def print_summary(df: pd.DataFrame):
    """
    Print a dataset summary.
    """

    valid = df["valid"].sum()
    invalid = len(df) - valid

    print("=" * 60)
    print("MatterGen Dataset Summary")
    print("=" * 60)
    print()

    print(f"Structures           : {len(df)}")
    print(f"Successfully loaded  : {valid}")
    print(f"Failed               : {invalid}")
    print()

    numeric_columns = [
        "num_atoms",
        "volume",
        "density",
        "a",
        "b",
        "c",
    ]

    for column in numeric_columns:

        values = df[column].dropna()

        if len(values) == 0:
            continue

        print(column)

        print(f"  Mean   : {values.mean():.3f}")
        print(f"  Median : {values.median():.3f}")
        print(f"  Std    : {values.std():.3f}")
        print(f"  Min    : {values.min():.3f}")
        print(f"  Max    : {values.max():.3f}")

        print()

    print("Most common space groups")

    counts = Counter(df["space_group"])

    for group, count in counts.most_common(10):

        print(f"{group:12s} {count}")

    print()


# -----------------------------------------------------------------------------


def save_csv(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Save the per-structure statistics.
    """

    csv_path = output_dir / "basic_statistics.csv"

    df.to_csv(
        csv_path,
        index=False,
    )

    summary = []

    numeric_columns = [
        "num_atoms",
        "volume",
        "density",
        "a",
        "b",
        "c",
        "alpha",
        "beta",
        "gamma",
    ]

    summary.append(
        {
            "Metric": "Structures",
            "Value": len(df),
        }
    )

    summary.append(
        {
            "Metric": "Valid structures",
            "Value": int(df["valid"].sum()),
        }
    )

    summary.append(
        {
            "Metric": "Invalid structures",
            "Value": int((~df["valid"]).sum()),
        }
    )

    for column in numeric_columns:

        values = df[column].dropna()

        if len(values) == 0:
            continue

        summary.extend([
            {
                "Metric": f"{column}_mean",
                "Value": values.mean(),
            },
            {
                "Metric": f"{column}_median",
                "Value": values.median(),
            },
            {
                "Metric": f"{column}_std",
                "Value": values.std(),
            },
        ])

    pd.DataFrame(summary).to_csv(
        output_dir / "summary.csv",
        index=False,
    )

    failed = df[df["valid"] == False]

    if not failed.empty:

        columns = ["file"]

        if "error" in failed.columns:
            columns.append("error")

        failed[columns].to_csv(
            output_dir / "failed_files.csv",
            index=False,
        )

    # -----------------------------------------------------------------------------


def plot_histogram(
    df: pd.DataFrame,
    column: str,
    xlabel: str,
    output_dir: Path,
):
    """
    Create histogram for a numeric column.
    """

    values = df[column].dropna()

    if len(values) == 0:
        return

    plt.figure()

    sns.histplot(
        values,
        bins=30,
        kde=True,
    )

    plt.title(xlabel)
    plt.xlabel(xlabel)
    plt.ylabel("Count")

    plt.tight_layout()

    plt.savefig(
        output_dir / f"{column}_histogram.pdf",
        bbox_inches="tight",
    )

    plt.savefig(
        output_dir / f"{column}_histogram.png",
        bbox_inches="tight",
    )

    plt.close()


# -----------------------------------------------------------------------------


def plot_spacegroups(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot most common space groups.
    """

    counts = (
        df["space_group"]
        .value_counts()
        .head(20)
    )

    plt.figure(figsize=(10, 6))

    sns.barplot(
        x=counts.values,
        y=counts.index,
    )

    plt.xlabel("Count")
    plt.ylabel("Space group")
    plt.title("Most common space groups")

    plt.tight_layout()

    plt.savefig(
        output_dir / "spacegroup_distribution.pdf",
        bbox_inches="tight",
    )

    plt.savefig(
        output_dir / "spacegroup_distribution.png",
        bbox_inches="tight",
    )

    plt.close()


# -----------------------------------------------------------------------------


def plot_formula_distribution(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot most common reduced formulas.
    """

    counts = (
        df["reduced_formula"]
        .value_counts()
        .head(20)
    )

    plt.figure(figsize=(10, 6))

    sns.barplot(
        x=counts.values,
        y=counts.index,
    )

    plt.xlabel("Count")
    plt.ylabel("Formula")
    plt.title("Most common compositions")

    plt.tight_layout()

    plt.savefig(
        output_dir / "composition_distribution.pdf",
        bbox_inches="tight",
    )

    plt.savefig(
        output_dir / "composition_distribution.png",
        bbox_inches="tight",
    )

    plt.close()


# -----------------------------------------------------------------------------


def plot_lattice_parameters(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot a, b and c distributions.
    """

    plt.figure(figsize=(10, 6))

    for column in ["a", "b", "c"]:

        sns.histplot(
            df[column].dropna(),
            bins=30,
            kde=True,
            label=column,
            stat="density",
            element="step",
            fill=False,
        )

    plt.xlabel("Lattice parameter (Å)")
    plt.ylabel("Density")
    plt.title("Lattice parameter distributions")

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        output_dir / "lattice_parameters.pdf",
        bbox_inches="tight",
    )

    plt.savefig(
        output_dir / "lattice_parameters.png",
        bbox_inches="tight",
    )

    plt.close()


# -----------------------------------------------------------------------------


def create_all_plots(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Generate all figures.
    """

    plot_histogram(
        df,
        "num_atoms",
        "Atoms per cell",
        output_dir,
    )

    plot_histogram(
        df,
        "volume",
        "Cell volume (Å³)",
        output_dir,
    )

    plot_histogram(
        df,
        "density",
        "Density (g/cm³)",
        output_dir,
    )

    plot_lattice_parameters(
        df,
        output_dir,
    )

    plot_spacegroups(
        df,
        output_dir,
    )

    plot_formula_distribution(
        df,
        output_dir,
    )


# -----------------------------------------------------------------------------


def main():

    args = parse_arguments()

    ensure_output_directory(
        args.output,
    )

    print()

    print("Searching for structures...")

    files = find_structure_files(
        args.input,
        args.recursive,
    )

    print(f"Found {len(files)} structures.")

    if len(files) == 0:

        print("No supported files found.")

        return

    print()

    df = analyze_dataset(
        files,
        args.workers,
    )

    print_summary(df)

    save_csv(
        df,
        args.output,
    )

    print("Generating plots...")

    create_all_plots(
        df,
        args.output,
    )

    print()

    print(f"Results written to {args.output.resolve()}")

    print("Done.")


# -----------------------------------------------------------------------------


if __name__ == "__main__":
    main()
