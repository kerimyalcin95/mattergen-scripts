#!/usr/bin/env python3
"""
02_symmetry_analysis.py

Analyze the crystallographic symmetry of a dataset of crystal structures.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes:
    * Space group number
    * Space group symbol
    * Crystal system
    * Bravais lattice
    * Point group
    * Hall symbol
    * Number of symmetry operations
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

DEFAULT_SYMPREC = 0.1

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
        description="Analyze crystallographic symmetry."
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

    parser.add_argument(
        "--symprec",
        default=DEFAULT_SYMPREC,
        type=float,
        help="Symmetry tolerance.",
    )

    return parser.parse_args()


# -----------------------------------------------------------------------------


def find_structure_files(
    folder: Path,
    recursive: bool,
):
    """
    Find supported structure files.
    """

    iterator = folder.rglob("*") if recursive else folder.glob("*")

    files = [
        path
        for path in iterator
        if path.suffix.lower() in SUPPORTED_EXTENSIONS
    ]

    files.sort()

    return files


# -----------------------------------------------------------------------------


def load_structure(
    path: Path,
):
    """
    Load a crystal structure.

    Parameters
    ----------
    path

    Returns
    -------
    pymatgen Structure
    """

    if path.suffix.lower() == ".cif":
        return Structure.from_file(path)

    atoms = ase_read(path)

    return AseAtomsAdaptor.get_structure(atoms)


# -----------------------------------------------------------------------------


def ensure_output_directory(
    output_dir: Path,
):
    """
    Create output directory.
    """

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


# -----------------------------------------------------------------------------


def create_analyzer(
    structure: Structure,
    symprec: float,
):
    """
    Create a SpacegroupAnalyzer.

    Returns
    -------
    SpacegroupAnalyzer
    """

    return SpacegroupAnalyzer(
        structure,
        symprec=symprec,
    )


# -----------------------------------------------------------------------------


def analyze_structure(
    path: Path,
    symprec: float,
):
    """
    Analyze symmetry of one structure.

    Parameters
    ----------
    path

    symprec

    Returns
    -------
    dict
    """

    try:

        structure = load_structure(path)

        analyzer = create_analyzer(
            structure,
            symprec,
        )

        dataset = analyzer.get_symmetry_dataset()

        result = {
            "file": path.name,
            "space_group_number": analyzer.get_space_group_number(),
            "space_group_symbol": analyzer.get_space_group_symbol(),
            "crystal_system": analyzer.get_crystal_system(),
            "bravais_lattice": analyzer.get_lattice_type(),
            "point_group": analyzer.get_point_group_symbol(),
            "hall_symbol": dataset.hall,
            "symmetry_operations": len(
                analyzer.get_symmetry_operations()
            ),
            "symprec": symprec,
            "valid": True,
        }

        return result

    except Exception as exc:

        return {
            "file": path.name,
            "space_group_number": np.nan,
            "space_group_symbol": "Invalid",
            "crystal_system": "Invalid",
            "bravais_lattice": "Invalid",
            "point_group": "Invalid",
            "hall_symbol": "Invalid",
            "symmetry_operations": np.nan,
            "symprec": symprec,
            "valid": False,
            "error": str(exc),
        }


# -----------------------------------------------------------------------------


def analyze_dataset(
    files,
    workers,
    symprec,
):
    """
    Analyze all structures.
    """

    iterator = tqdm(
        files,
        desc="Analyzing symmetry",
    )

    rows = Parallel(
        n_jobs=workers,
    )(
        delayed(analyze_structure)(
            path,
            symprec,
        )
        for path in iterator
    )

    return pd.DataFrame(rows)

# -----------------------------------------------------------------------------


def print_summary(
    df: pd.DataFrame,
):
    """
    Print a dataset summary.
    """

    valid = int(df["valid"].sum())
    invalid = len(df) - valid

    print("=" * 60)
    print("Symmetry Analysis Summary")
    print("=" * 60)
    print()

    print(f"Structures           : {len(df)}")
    print(f"Successfully analyzed: {valid}")
    print(f"Failed               : {invalid}")
    print()

    if valid == 0:
        return

    print(f"Unique space groups  : {df['space_group_symbol'].nunique()}")
    print(f"Unique point groups  : {df['point_group'].nunique()}")
    print(f"Unique Hall symbols  : {df['hall_symbol'].nunique()}")
    print()

    print("Crystal systems")

    counts = (
        df.loc[df["valid"], "crystal_system"]
        .value_counts()
        .sort_index()
    )

    total = counts.sum()

    for system, count in counts.items():

        percentage = 100 * count / total

        print(f"{system:15s} {count:8d} ({percentage:5.1f}%)")

    print()

    print("Top 10 space groups")

    counts = Counter(
        df.loc[df["valid"], "space_group_symbol"]
    )

    for group, count in counts.most_common(10):

        print(f"{group:12s} {count}")

    print()


# -----------------------------------------------------------------------------


def save_csv(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Save analysis results.
    """

    df.to_csv(
        output_dir / "symmetry_analysis.csv",
        index=False,
    )

    summary = []

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

    valid = df[df["valid"]]

    if not valid.empty:

        summary.extend(
            [
                {
                    "Metric": "Unique space groups",
                    "Value": valid["space_group_symbol"].nunique(),
                },
                {
                    "Metric": "Unique point groups",
                    "Value": valid["point_group"].nunique(),
                },
                {
                    "Metric": "Unique Hall symbols",
                    "Value": valid["hall_symbol"].nunique(),
                },
                {
                    "Metric": "Unique crystal systems",
                    "Value": valid["crystal_system"].nunique(),
                },
                {
                    "Metric": "Mean symmetry operations",
                    "Value": valid["symmetry_operations"].mean(),
                },
                {
                    "Metric": "Median symmetry operations",
                    "Value": valid["symmetry_operations"].median(),
                },
                {
                    "Metric": "Maximum symmetry operations",
                    "Value": valid["symmetry_operations"].max(),
                },
            ]
        )

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


def create_output_directory(
    base_output: Path,
):
    """
    Create the dedicated symmetry output directory.

    Returns
    -------
    Path
    """

    output_dir = base_output / "symmetry"

    ensure_output_directory(
        output_dir,
    )

    return output_dir
# -----------------------------------------------------------------------------


def plot_crystal_systems(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot crystal system distribution.
    """

    counts = (
        df.loc[df["valid"], "crystal_system"]
        .value_counts()
        .sort_values(ascending=False)
    )

    if counts.empty:
        return

    plt.figure(figsize=(8, 5))

    sns.barplot(
        x=counts.values,
        y=counts.index,
    )

    plt.xlabel("Count")
    plt.ylabel("Crystal system")
    plt.title("Crystal system distribution")

    plt.tight_layout()

    for extension in ("pdf", "png"):
        plt.savefig(
            output_dir / f"crystal_system_distribution.{extension}",
            bbox_inches="tight",
        )

    plt.close()


# -----------------------------------------------------------------------------


def plot_space_groups(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot the 20 most common space groups.
    """

    counts = (
        df.loc[df["valid"], "space_group_symbol"]
        .value_counts()
        .head(20)
    )

    if counts.empty:
        return

    plt.figure(figsize=(10, 6))

    sns.barplot(
        x=counts.values,
        y=counts.index,
    )

    plt.xlabel("Count")
    plt.ylabel("Space group")
    plt.title("Most common space groups")

    plt.tight_layout()

    for extension in ("pdf", "png"):
        plt.savefig(
            output_dir / f"space_group_distribution.{extension}",
            bbox_inches="tight",
        )

    plt.close()


# -----------------------------------------------------------------------------


def plot_space_group_numbers(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot histogram of space group numbers.
    """

    values = df.loc[
        df["valid"],
        "space_group_number",
    ].dropna()

    if values.empty:
        return

    plt.figure()

    sns.histplot(
        values,
        bins=230,
        kde=False,
    )

    plt.xlabel("Space group number")
    plt.ylabel("Count")
    plt.title("Space group number distribution")

    plt.tight_layout()

    for extension in ("pdf", "png"):
        plt.savefig(
            output_dir / f"space_group_number_histogram.{extension}",
            bbox_inches="tight",
        )

    plt.close()


# -----------------------------------------------------------------------------


def plot_point_groups(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot point group distribution.
    """

    counts = (
        df.loc[df["valid"], "point_group"]
        .value_counts()
    )

    if counts.empty:
        return

    plt.figure(figsize=(10, 6))

    sns.barplot(
        x=counts.values,
        y=counts.index,
    )

    plt.xlabel("Count")
    plt.ylabel("Point group")
    plt.title("Point group distribution")

    plt.tight_layout()

    for extension in ("pdf", "png"):
        plt.savefig(
            output_dir / f"point_group_distribution.{extension}",
            bbox_inches="tight",
        )

    plt.close()


# -----------------------------------------------------------------------------


def plot_symmetry_operations(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot symmetry operation distribution.
    """

    values = df.loc[
        df["valid"],
        "symmetry_operations",
    ].dropna()

    if values.empty:
        return

    plt.figure()

    sns.histplot(
        values,
        bins=30,
        kde=True,
    )

    plt.xlabel("Number of symmetry operations")
    plt.ylabel("Count")
    plt.title("Symmetry operations")

    plt.tight_layout()

    for extension in ("pdf", "png"):
        plt.savefig(
            output_dir / f"symmetry_operations_histogram.{extension}",
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

    plot_crystal_systems(
        df,
        output_dir,
    )

    plot_space_groups(
        df,
        output_dir,
    )

    plot_space_group_numbers(
        df,
        output_dir,
    )

    plot_point_groups(
        df,
        output_dir,
    )

    plot_symmetry_operations(
        df,
        output_dir,
    )

# -----------------------------------------------------------------------------

def main():

    args = parse_arguments()

    output_dir = create_output_directory(
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
        files=files,
        workers=args.workers,
        symprec=args.symprec,
    )

    print_summary(df)

    save_csv(
        df,
        output_dir,
    )

    print("Generating plots...")

    create_all_plots(
        df,
        output_dir,
    )

    print()
    print(f"Results written to {output_dir.resolve()}")
    print("Done.")


# -----------------------------------------------------------------------------


if __name__ == "__main__":
    main()
