#!/usr/bin/env python3
"""
03_coordination_analysis.py

Analyze local coordination environments of crystal structures using both
CrystalNN and VoronoiNN.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes per-atom coordination numbers with:
    * CrystalNN
    * VoronoiNN
- Computes per-structure statistics:
    * mean, median, std, min, max coordination
    * agreement fraction between both methods
    * mean and max absolute difference
- Saves:
    * CSV with one row per structure
    * Summary CSV
    * Failed files CSV
    * Publication-quality PDF/PNG figures
- Prints a human-readable summary to terminal

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
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from joblib import Parallel, delayed
from tqdm import tqdm

import matplotlib.pyplot as plt
import seaborn as sns

from ase.io import read as ase_read
from pymatgen.analysis.local_env import CrystalNN, VoronoiNN
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

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
        description="Analyze coordination environments using CrystalNN and VoronoiNN."
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
        help="Base output directory.",
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


def find_structure_files(folder: Path, recursive: bool):
    """
    Find supported structure files.
    """

    iterator = folder.rglob("*") if recursive else folder.glob("*")

    files = [path for path in iterator if path.suffix.lower() in SUPPORTED_EXTENSIONS]
    files.sort()
    return files


# -----------------------------------------------------------------------------


def load_structure(path: Path):
    """
    Load a crystal structure using pymatgen or ASE.
    """

    if path.suffix.lower() == ".cif":
        return Structure.from_file(path)

    atoms = ase_read(path)
    return AseAtomsAdaptor.get_structure(atoms)


# -----------------------------------------------------------------------------


def ensure_output_directory(output_dir: Path):
    """
    Create output directory if necessary.
    """

    output_dir.mkdir(parents=True, exist_ok=True)


# -----------------------------------------------------------------------------


def create_output_directory(base_output: Path):
    """
    Create the dedicated coordination output directory.
    """

    output_dir = base_output / "coordination"
    ensure_output_directory(output_dir)
    return output_dir


# -----------------------------------------------------------------------------


def coordination_number_from_nninfo(nn_info):
    """
    Convert neighbor information into an integer coordination number.
    """

    if nn_info is None:
        return 0
    return len(nn_info)


# -----------------------------------------------------------------------------


def analyze_structure(path: Path):
    """
    Analyze coordination for one structure.
    """

    try:
        structure = load_structure(path)
        crystal_nn = crystal_nn = CrystalNN(
            distance_cutoffs=None,
            x_diff_weight=0,
            porous_adjustment=False,
        )
        voronoi_nn = VoronoiNN()

        crystal_cns = []
        voronoi_cns = []
        per_site_differences = []
        agreement_hits = 0

        for site_index in range(len(structure)):
            try:
                c_nn = crystal_nn.get_nn_info(structure, site_index)
                c_cn = coordination_number_from_nninfo(c_nn)
            except Exception:
                c_cn = np.nan

            try:
                v_nn = voronoi_nn.get_nn_info(structure, site_index)
                v_cn = coordination_number_from_nninfo(v_nn)
            except Exception:
                v_cn = np.nan

            crystal_cns.append(c_cn)
            voronoi_cns.append(v_cn)

            if np.isfinite(c_cn) and np.isfinite(v_cn):
                diff = abs(float(c_cn) - float(v_cn))
                per_site_differences.append(diff)
                if diff == 0:
                    agreement_hits += 1

        crystal_values = pd.Series(crystal_cns, dtype="float64").dropna()
        voronoi_values = pd.Series(voronoi_cns, dtype="float64").dropna()
        diff_values = pd.Series(per_site_differences, dtype="float64").dropna()

        formula = structure.composition.formula
        reduced_formula = structure.composition.reduced_formula

        result = {
            "file": path.name,
            "formula": formula,
            "reduced_formula": reduced_formula,
            "num_atoms": len(structure),
            "crystalnn_mean_coordination": crystal_values.mean() if not crystal_values.empty else np.nan,
            "crystalnn_median_coordination": crystal_values.median() if not crystal_values.empty else np.nan,
            "crystalnn_std_coordination": crystal_values.std() if not crystal_values.empty else np.nan,
            "crystalnn_min_coordination": crystal_values.min() if not crystal_values.empty else np.nan,
            "crystalnn_max_coordination": crystal_values.max() if not crystal_values.empty else np.nan,
            "voronoi_mean_coordination": voronoi_values.mean() if not voronoi_values.empty else np.nan,
            "voronoi_median_coordination": voronoi_values.median() if not voronoi_values.empty else np.nan,
            "voronoi_std_coordination": voronoi_values.std() if not voronoi_values.empty else np.nan,
            "voronoi_min_coordination": voronoi_values.min() if not voronoi_values.empty else np.nan,
            "voronoi_max_coordination": voronoi_values.max() if not voronoi_values.empty else np.nan,
            "mean_absolute_difference": diff_values.mean() if not diff_values.empty else np.nan,
            "max_absolute_difference": diff_values.max() if not diff_values.empty else np.nan,
            "agreement_fraction": agreement_hits / len(structure) if len(structure) > 0 else np.nan,
            "valid": True,
        }

        return result

    except Exception as exc:
        return {
            "file": path.name,
            "formula": None,
            "reduced_formula": None,
            "num_atoms": np.nan,
            "crystalnn_mean_coordination": np.nan,
            "crystalnn_median_coordination": np.nan,
            "crystalnn_std_coordination": np.nan,
            "crystalnn_min_coordination": np.nan,
            "crystalnn_max_coordination": np.nan,
            "voronoi_mean_coordination": np.nan,
            "voronoi_median_coordination": np.nan,
            "voronoi_std_coordination": np.nan,
            "voronoi_min_coordination": np.nan,
            "voronoi_max_coordination": np.nan,
            "mean_absolute_difference": np.nan,
            "max_absolute_difference": np.nan,
            "agreement_fraction": np.nan,
            "valid": False,
            "error": str(exc),
        }


# -----------------------------------------------------------------------------


def analyze_dataset(files, workers):
    """
    Analyze all structures in parallel.
    """

    iterator = tqdm(files, desc="Analyzing coordination")

    rows = Parallel(n_jobs=workers)(
        delayed(analyze_structure)(path) for path in iterator
    )

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------


def print_summary(df: pd.DataFrame):
    """
    Print a dataset summary.
    """

    valid = int(df["valid"].sum())
    invalid = len(df) - valid

    print("=" * 60)
    print("Coordination Analysis Summary")
    print("=" * 60)
    print()

    print(f"Structures           : {len(df)}")
    print(f"Successfully analyzed : {valid}")
    print(f"Failed               : {invalid}")
    print()

    if valid == 0:
        return

    metrics = [
        "crystalnn_mean_coordination",
        "voronoi_mean_coordination",
        "mean_absolute_difference",
        "agreement_fraction",
    ]

    for column in metrics:
        values = df[column].dropna()
        if values.empty:
            continue
        print(column)
        print(f"  Mean   : {values.mean():.3f}")
        print(f"  Median : {values.median():.3f}")
        print(f"  Std    : {values.std():.3f}")
        print(f"  Min    : {values.min():.3f}")
        print(f"  Max    : {values.max():.3f}")
        print()

    print("Most common formulas")
    counts = Counter(df.loc[df["valid"], "reduced_formula"])
    for formula, count in counts.most_common(10):
        print(f"{formula:20s} {count}")
    print()


# -----------------------------------------------------------------------------


def save_csv(df: pd.DataFrame, output_dir: Path):
    """
    Save analysis results.
    """

    df.to_csv(output_dir / "coordination_analysis.csv", index=False)

    summary = [
        {"Metric": "Structures", "Value": len(df)},
        {"Metric": "Valid structures", "Value": int(df["valid"].sum())},
        {"Metric": "Invalid structures", "Value": int((~df["valid"]).sum())},
    ]

    valid = df[df["valid"]]
    numeric_columns = [
        "crystalnn_mean_coordination",
        "voronoi_mean_coordination",
        "mean_absolute_difference",
        "agreement_fraction",
    ]

    for column in numeric_columns:
        values = valid[column].dropna()
        if values.empty:
            continue
        summary.extend([
            {"Metric": f"{column}_mean", "Value": values.mean()},
            {"Metric": f"{column}_median", "Value": values.median()},
            {"Metric": f"{column}_std", "Value": values.std()},
        ])

    pd.DataFrame(summary).to_csv(output_dir / "summary.csv", index=False)

    failed = df[df["valid"] == False]
    if not failed.empty:
        columns = ["file"]
        if "error" in failed.columns:
            columns.append("error")
        failed[columns].to_csv(output_dir / "failed_files.csv", index=False)


# -----------------------------------------------------------------------------


def plot_histogram(df: pd.DataFrame, column: str, xlabel: str, output_dir: Path, filename: str):
    """
    Create histogram for a numeric column.
    """

    values = df[column].dropna()
    if values.empty:
        return

    plt.figure()
    sns.histplot(values, bins=30, kde=True)
    plt.title(xlabel)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.tight_layout()

    for ext in ("pdf", "png"):
        plt.savefig(output_dir / f"{filename}.{ext}", bbox_inches="tight")
    plt.close()


# -----------------------------------------------------------------------------


def plot_comparison(df: pd.DataFrame, output_dir: Path):
    """
    Scatter plot comparing CrystalNN and VoronoiNN.
    """

    valid = df[df["valid"]].dropna(subset=["crystalnn_mean_coordination", "voronoi_mean_coordination"])
    if valid.empty:
        return

    plt.figure(figsize=(6, 6))
    sns.scatterplot(
        data=valid,
        x="crystalnn_mean_coordination",
        y="voronoi_mean_coordination",
        s=40,
    )
    max_value = max(
        valid["crystalnn_mean_coordination"].max(),
        valid["voronoi_mean_coordination"].max(),
    )
    plt.plot([0, max_value], [0, max_value], linestyle="--")
    plt.xlabel("CrystalNN mean coordination")
    plt.ylabel("VoronoiNN mean coordination")
    plt.title("CrystalNN vs VoronoiNN")
    plt.tight_layout()

    for ext in ("pdf", "png"):
        plt.savefig(output_dir / f"coordination_scatter.{ext}", bbox_inches="tight")
    plt.close()


# -----------------------------------------------------------------------------


def plot_difference(df: pd.DataFrame, output_dir: Path):
    """
    Plot distribution of coordination differences.
    """

    values = df.loc[df["valid"], "mean_absolute_difference"].dropna()
    if values.empty:
        return

    plt.figure()
    sns.histplot(values, bins=30, kde=True)
    plt.xlabel("Mean absolute difference")
    plt.ylabel("Count")
    plt.title("CrystalNN vs VoronoiNN difference")
    plt.tight_layout()

    for ext in ("pdf", "png"):
        plt.savefig(output_dir / f"coordination_difference_histogram.{ext}", bbox_inches="tight")
    plt.close()


# -----------------------------------------------------------------------------


def plot_boxplot(df: pd.DataFrame, output_dir: Path):
    """
    Compare coordination number distributions.
    """

    valid = df[df["valid"]]
    if valid.empty:
        return

    melted = pd.melt(
        valid,
        value_vars=["crystalnn_mean_coordination", "voronoi_mean_coordination"],
        var_name="method",
        value_name="coordination",
    )

    plt.figure(figsize=(7, 5))
    sns.boxplot(data=melted, x="method", y="coordination")
    plt.xlabel("Method")
    plt.ylabel("Mean coordination")
    plt.title("Coordination comparison")
    plt.tight_layout()

    for ext in ("pdf", "png"):
        plt.savefig(output_dir / f"coordination_boxplot.{ext}", bbox_inches="tight")
    plt.close()


# -----------------------------------------------------------------------------


def plot_formula_distribution(df: pd.DataFrame, output_dir: Path):
    """
    Plot most common reduced formulas.
    """

    counts = df.loc[df["valid"], "reduced_formula"].value_counts().head(20)
    if counts.empty:
        return

    plt.figure(figsize=(10, 6))
    sns.barplot(x=counts.values, y=counts.index)
    plt.xlabel("Count")
    plt.ylabel("Formula")
    plt.title("Most common compositions")
    plt.tight_layout()

    for ext in ("pdf", "png"):
        plt.savefig(output_dir / f"composition_distribution.{ext}", bbox_inches="tight")
    plt.close()


# -----------------------------------------------------------------------------


def create_all_plots(df: pd.DataFrame, output_dir: Path):
    """
    Generate all figures.
    """

    plot_histogram(df, "crystalnn_mean_coordination", "CrystalNN mean coordination", output_dir, "crystalnn_coordination_histogram")
    plot_histogram(df, "voronoi_mean_coordination", "VoronoiNN mean coordination", output_dir, "voronoi_coordination_histogram")
    plot_histogram(df, "agreement_fraction", "Agreement fraction", output_dir, "agreement_fraction_histogram")
    plot_difference(df, output_dir)
    plot_comparison(df, output_dir)
    plot_boxplot(df, output_dir)
    plot_formula_distribution(df, output_dir)


# -----------------------------------------------------------------------------


def main():
    args = parse_arguments()

    output_dir = create_output_directory(args.output)

    print()
    print("Searching for structures...")

    files = find_structure_files(args.input, args.recursive)
    print(f"Found {len(files)} structures.")

    if len(files) == 0:
        print("No supported files found.")
        return

    print()
    df = analyze_dataset(files, args.workers)

    print_summary(df)

    save_csv(df, output_dir)

    print("Generating plots...")
    create_all_plots(df, output_dir)

    print()
    print(f"Results written to {output_dir.resolve()}")
    print("Done.")


# -----------------------------------------------------------------------------


if __name__ == "__main__":
    main()
