#!/usr/bin/env python3
"""
04_bond_analysis.py

Analyze bond lengths, bond angles and nearest-neighbor statistics.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes:
    * Bond length statistics
    * Bond angle statistics
    * Nearest-neighbor statistics
    * Mean coordination number
- Saves:
    * CSV with one row per structure
    * Summary CSV
    * Failed files CSV
    * Publication-quality PDF and PNG figures
- Prints summary to terminal

Supported libraries:
    pymatgen
    ase
    numpy
    pandas
    matplotlib
    seaborn
    scipy
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

from scipy.spatial.distance import pdist

from ase.io import read as ase_read

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.local_env import CrystalNN

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
        files,
        args.workers,
    )

    print_summary(
        df,
    )

    save_csv(
        df,
        output_dir,
    )


    print("Creating figures...")

    create_all_plots(
        df,
        output_dir,
    )

    print()
    print(f"Results written to {output_dir.resolve()}")
    print("Done.")


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description="Analyze bond lengths, bond angles and nearest-neighbor statistics."
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
    Load a crystal structure using pymatgen or ASE.
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
    Create output directory if necessary.
    """

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


# -----------------------------------------------------------------------------


def create_output_directory(
    base_output: Path,
):
    """
    Create the dedicated bond output directory.

    Returns
    -------
    Path
    """

    output_dir = base_output / "bond"

    ensure_output_directory(
        output_dir,
    )

    return output_dir


# -----------------------------------------------------------------------------


def safe_series_statistics(
    values,
):
    """
    Compute common statistics for a numeric sequence.

    Returns
    -------
    tuple
        mean, median, std, min, max
    """

    values = pd.Series(
        values,
        dtype="float64",
    ).dropna()

    if values.empty:
        return (
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        )

    return (
        values.mean(),
        values.median(),
        values.std(),
        values.min(),
        values.max(),
    )


# -----------------------------------------------------------------------------


def unique_pair(
    index_a: int,
    index_b: int,
):
    """
    Return an order-independent atom pair.
    """

    return tuple(
        sorted(
            (
                index_a,
                index_b,
            )
        )
    )

# -----------------------------------------------------------------------------


def extract_bonds(
    structure: Structure,
):
    """
    Extract unique bonds using CrystalNN.

    Returns
    -------
    tuple[list[float], list[int], dict[int, list[int]]]
        bond_lengths
        coordination_numbers
        neighbor_map
    """

    crystal_nn = CrystalNN(
        distance_cutoffs=None,
        x_diff_weight=0,
        porous_adjustment=False,
    )

    unique_bonds = set()

    bond_lengths = []

    coordination_numbers = []

    neighbor_map = {
        index: []
        for index in range(len(structure))
    }

    for site_index in range(len(structure)):

        try:

            neighbors = crystal_nn.get_nn_info(
                structure,
                site_index,
            )

        except Exception:

            coordination_numbers.append(np.nan)

            continue

        coordination_numbers.append(
            len(neighbors)
        )

        for neighbor in neighbors:

            neighbor_index = neighbor["site_index"]

            pair = unique_pair(
                site_index,
                neighbor_index,
            )

            neighbor_map[site_index].append(
                neighbor_index
            )

            if pair in unique_bonds:
                continue

            unique_bonds.add(pair)

            distance = structure.get_distance(
                site_index,
                neighbor_index,
            )

            bond_lengths.append(
                float(distance)
            )

    return (
        bond_lengths,
        coordination_numbers,
        neighbor_map,
    )


# -----------------------------------------------------------------------------


def compute_nearest_neighbor_distances(
    structure: Structure,
    neighbor_map,
):
    """
    Compute nearest-neighbor distance for every atom.

    Returns
    -------
    list[float]
    """

    nearest_distances = []

    for atom_index in range(len(structure)):

        neighbors = neighbor_map.get(
            atom_index,
            [],
        )

        if len(neighbors) == 0:

            nearest_distances.append(
                np.nan,
            )

            continue

        distances = [
            structure.get_distance(
                atom_index,
                neighbor_index,
            )
            for neighbor_index in neighbors
        ]

        nearest_distances.append(
            float(min(distances))
        )

    return nearest_distances


# -----------------------------------------------------------------------------


def summarize_bond_statistics(
    bond_lengths,
    nearest_neighbor_distances,
    coordination_numbers,
):
    """
    Compute summary statistics.

    Returns
    -------
    dict
    """

    (
        bond_mean,
        bond_median,
        bond_std,
        bond_min,
        bond_max,
    ) = safe_series_statistics(
        bond_lengths,
    )

    (
        nn_mean,
        nn_median,
        nn_std,
        nn_min,
        nn_max,
    ) = safe_series_statistics(
        nearest_neighbor_distances,
    )

    (
        coordination_mean,
        coordination_median,
        coordination_std,
        coordination_min,
        coordination_max,
    ) = safe_series_statistics(
        coordination_numbers,
    )

    return {
        "bond_mean": bond_mean,
        "bond_median": bond_median,
        "bond_std": bond_std,
        "bond_min": bond_min,
        "bond_max": bond_max,
        "nearest_neighbor_mean": nn_mean,
        "nearest_neighbor_median": nn_median,
        "nearest_neighbor_std": nn_std,
        "nearest_neighbor_min": nn_min,
        "nearest_neighbor_max": nn_max,
        "coordination_mean": coordination_mean,
        "coordination_median": coordination_median,
        "coordination_std": coordination_std,
        "coordination_min": coordination_min,
        "coordination_max": coordination_max,
        "number_of_bonds": len(bond_lengths),
    }
# -----------------------------------------------------------------------------


def compute_bond_angles(
    structure: Structure,
    neighbor_map,
):
    """
    Compute all bond angles in the structure.

    For every central atom B and neighboring atoms
    A and C, the angle A-B-C is calculated.

    Returns
    -------
    list[float]
        Bond angles in degrees.
    """

    bond_angles = []

    for center_index in range(len(structure)):

        neighbors = neighbor_map.get(
            center_index,
            [],
        )

        if len(neighbors) < 2:
            continue

        center = structure[center_index].coords

        for i in range(len(neighbors) - 1):

            first_index = neighbors[i]

            vector_1 = (
                structure[first_index].coords
                - center
            )

            norm_1 = np.linalg.norm(
                vector_1,
            )

            if norm_1 == 0:
                continue

            for j in range(i + 1, len(neighbors)):

                second_index = neighbors[j]

                vector_2 = (
                    structure[second_index].coords
                    - center
                )

                norm_2 = np.linalg.norm(
                    vector_2,
                )

                if norm_2 == 0:
                    continue

                cosine = np.dot(
                    vector_1,
                    vector_2,
                ) / (
                    norm_1 * norm_2
                )

                cosine = np.clip(
                    cosine,
                    -1.0,
                    1.0,
                )

                angle = np.degrees(
                    np.arccos(
                        cosine,
                    )
                )

                bond_angles.append(
                    float(angle)
                )

    return bond_angles


# -----------------------------------------------------------------------------


def summarize_angle_statistics(
    bond_angles,
):
    """
    Compute bond-angle statistics.

    Returns
    -------
    dict
    """

    (
        angle_mean,
        angle_median,
        angle_std,
        angle_min,
        angle_max,
    ) = safe_series_statistics(
        bond_angles,
    )

    return {
        "bond_angle_mean": angle_mean,
        "bond_angle_median": angle_median,
        "bond_angle_std": angle_std,
        "bond_angle_min": angle_min,
        "bond_angle_max": angle_max,
        "number_of_angles": len(
            bond_angles,
        ),
    }


# -----------------------------------------------------------------------------


def collect_bond_information(
    structure: Structure,
):
    """
    Compute all bond-related information for one structure.

    Returns
    -------
    dict
    """

    (
        bond_lengths,
        coordination_numbers,
        neighbor_map,
    ) = extract_bonds(
        structure,
    )

    nearest_neighbor_distances = (
        compute_nearest_neighbor_distances(
            structure,
            neighbor_map,
        )
    )

    bond_angles = compute_bond_angles(
        structure,
        neighbor_map,
    )

    results = summarize_bond_statistics(
        bond_lengths,
        nearest_neighbor_distances,
        coordination_numbers,
    )

    results.update(
        summarize_angle_statistics(
            bond_angles,
        )
    )

    return results
# -----------------------------------------------------------------------------


def analyze_structure(
    path: Path,
):
    """
    Analyze one crystal structure.

    Parameters
    ----------
    path

    Returns
    -------
    dict
    """

    try:

        structure = load_structure(
            path,
        )

        statistics = collect_bond_information(
            structure,
        )

        result = {
            "file": path.name,
            "formula": structure.composition.formula,
            "reduced_formula": structure.composition.reduced_formula,
            "num_atoms": len(structure),
            **statistics,
            "valid": True,
        }

        return result

    except Exception as exc:

        return {
            "file": path.name,
            "formula": None,
            "reduced_formula": None,
            "num_atoms": np.nan,
            "bond_mean": np.nan,
            "bond_median": np.nan,
            "bond_std": np.nan,
            "bond_min": np.nan,
            "bond_max": np.nan,
            "nearest_neighbor_mean": np.nan,
            "nearest_neighbor_median": np.nan,
            "nearest_neighbor_std": np.nan,
            "nearest_neighbor_min": np.nan,
            "nearest_neighbor_max": np.nan,
            "coordination_mean": np.nan,
            "coordination_median": np.nan,
            "coordination_std": np.nan,
            "coordination_min": np.nan,
            "coordination_max": np.nan,
            "bond_angle_mean": np.nan,
            "bond_angle_median": np.nan,
            "bond_angle_std": np.nan,
            "bond_angle_min": np.nan,
            "bond_angle_max": np.nan,
            "number_of_bonds": np.nan,
            "number_of_angles": np.nan,
            "valid": False,
            "error": str(exc),
        }


# -----------------------------------------------------------------------------


def analyze_dataset(
    files,
    workers,
):
    """
    Analyze an entire dataset.

    Parameters
    ----------
    files

    workers

    Returns
    -------
    pandas.DataFrame
    """

    iterator = tqdm(
        files,
        desc="Analyzing bonds",
    )

    rows = Parallel(
        n_jobs=workers,
    )(
        delayed(
            analyze_structure,
        )(path)
        for path in iterator
    )

    return pd.DataFrame(
        rows,
    )
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
    print("Bond Analysis Summary")
    print("=" * 60)
    print()

    print(f"Structures            : {len(df)}")
    print(f"Successfully analyzed : {valid}")
    print(f"Failed                : {invalid}")
    print()

    if valid == 0:
        return

    metrics = [
        "bond_mean",
        "bond_angle_mean",
        "nearest_neighbor_mean",
        "coordination_mean",
    ]

    for column in metrics:

        values = df.loc[
            df["valid"],
            column,
        ].dropna()

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

    counts = Counter(
        df.loc[
            df["valid"],
            "reduced_formula",
        ]
    )

    for formula, count in counts.most_common(10):

        print(f"{formula:20s} {count}")

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
        output_dir / "bond_analysis.csv",
        index=False,
    )

    summary = [
        {
            "Metric": "Structures",
            "Value": len(df),
        },
        {
            "Metric": "Valid structures",
            "Value": int(df["valid"].sum()),
        },
        {
            "Metric": "Invalid structures",
            "Value": int((~df["valid"]).sum()),
        },
    ]

    valid = df[
        df["valid"]
    ]

    numeric_columns = [
        "bond_mean",
        "bond_angle_mean",
        "nearest_neighbor_mean",
        "coordination_mean",
        "number_of_bonds",
        "number_of_angles",
    ]

    for column in numeric_columns:

        values = valid[
            column
        ].dropna()

        if values.empty:
            continue

        summary.extend(
            [
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
                {
                    "Metric": f"{column}_min",
                    "Value": values.min(),
                },
                {
                    "Metric": f"{column}_max",
                    "Value": values.max(),
                },
            ]
        )

    pd.DataFrame(
        summary,
    ).to_csv(
        output_dir / "summary.csv",
        index=False,
    )

    failed = df[
        df["valid"] == False
    ]

    if not failed.empty:

        columns = [
            "file",
        ]

        if "error" in failed.columns:
            columns.append(
                "error",
            )

        failed[
            columns
        ].to_csv(
            output_dir / "failed_files.csv",
            index=False,
        )

# -----------------------------------------------------------------------------


def save_figure(
    figure,
    output_dir: Path,
    filename: str,
):
    """
    Save a figure as PNG and PDF.
    """

    figure.tight_layout()

    figure.savefig(
        output_dir / f"{filename}.png",
        dpi=300,
    )

    figure.savefig(
        output_dir / f"{filename}.pdf",
    )

    plt.close(
        figure,
    )


# -----------------------------------------------------------------------------


def plot_histogram(
    values,
    xlabel: str,
    title: str,
    filename: str,
    output_dir: Path,
    bins: int = 50,
):
    """
    Plot a histogram.
    """

    values = pd.Series(values).dropna()

    if values.empty:
        return

    fig, ax = plt.subplots()

    sns.histplot(
        values,
        bins=bins,
        kde=True,
        stat="count",
        ax=ax,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_title(title)

    save_figure(
        fig,
        output_dir,
        filename,
    )


# -----------------------------------------------------------------------------


def plot_boxplots(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot boxplots of the main statistics.
    """

    columns = [
        "bond_mean",
        "bond_angle_mean",
        "nearest_neighbor_mean",
        "coordination_mean",
    ]

    labels = [
        "Bond\nLength",
        "Bond\nAngle",
        "Nearest\nNeighbor",
        "Coordination",
    ]

    plot_df = df[
        columns
    ].copy()

    if plot_df.dropna(how="all").empty:
        return

    fig, ax = plt.subplots(
        figsize=(8, 5),
    )

    sns.boxplot(
        data=plot_df,
        ax=ax,
    )

    ax.set_xticklabels(
        labels,
    )

    ax.set_ylabel("Value")

    ax.set_title(
        "Bond Statistics",
    )

    save_figure(
        fig,
        output_dir,
        "bond_statistics_boxplot",
    )


# -----------------------------------------------------------------------------


def create_all_plots(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Create every figure.
    """

    valid = df[
        df["valid"]
    ]

    if valid.empty:
        return

    plot_histogram(
        valid["bond_mean"],
        xlabel="Mean Bond Length (Å)",
        title="Bond Length Distribution",
        filename="bond_length_histogram",
        output_dir=output_dir,
    )

    plot_histogram(
        valid["bond_angle_mean"],
        xlabel="Mean Bond Angle (°)",
        title="Bond Angle Distribution",
        filename="bond_angle_histogram",
        output_dir=output_dir,
    )

    plot_histogram(
        valid["nearest_neighbor_mean"],
        xlabel="Nearest Neighbor Distance (Å)",
        title="Nearest Neighbor Distribution",
        filename="nearest_neighbor_histogram",
        output_dir=output_dir,
    )

    plot_histogram(
        valid["coordination_mean"],
        xlabel="Coordination Number",
        title="Coordination Number Distribution",
        filename="coordination_histogram",
        output_dir=output_dir,
        bins=20,
    )

    plot_boxplots(
        valid,
        output_dir,
    )


if __name__ == "__main__":
    main()
