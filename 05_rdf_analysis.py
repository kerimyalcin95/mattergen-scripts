#!/usr/bin/env python3
"""
05_rdf_analysis.py

Analyze radial distribution functions (RDF) of crystal structures.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes radial distribution functions
- Computes:
    * First RDF peak position
    * First RDF peak height
    * RDF integral
    * Maximum RDF
- Saves:
    * CSV with one row per structure
    * Average RDF CSV
    * Summary CSV
    * Failed files CSV
    * Publication-quality PDF and PNG figures
- Prints a dataset summary

Supported libraries
-------------------
- pymatgen
- ase
- numpy
- pandas
- scipy
- matplotlib
- seaborn
- tqdm
- joblib
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

from scipy.integrate import trapezoid
from scipy.signal import find_peaks

from ase.io import read as ase_read

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from pymatviz.rdf.helpers import calculate_rdf

# -----------------------------------------------------------------------------

SUPPORTED_EXTENSIONS = {
    ".cif",
    ".xyz",
    ".extxyz",
}

DEFAULT_R_MAX = 10.0
DEFAULT_BIN_SIZE = 0.02

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
        description="Analyze radial distribution functions."
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

    parser.add_argument(
        "--rmax",
        default=DEFAULT_R_MAX,
        type=float,
        help="Maximum RDF distance in Å.",
    )

    parser.add_argument(
        "--bin-size",
        default=DEFAULT_BIN_SIZE,
        type=float,
        help="Histogram bin width in Å.",
    )

    return parser.parse_args()


# -----------------------------------------------------------------------------


def find_structure_files(
    folder: Path,
    recursive: bool,
):
    """
    Find supported crystal structures.
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
    Load a structure.

    Returns
    -------
    pymatgen Structure
    """

    if path.suffix.lower() == ".cif":
        return Structure.from_file(path)

    atoms = ase_read(path)

    return AseAtomsAdaptor.get_structure(
        atoms,
    )


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


def create_output_directory(
    base_output: Path,
):
    """
    Create rdf output directory.
    """

    output_dir = base_output / "rdf"

    ensure_output_directory(
        output_dir,
    )

    return output_dir


# -----------------------------------------------------------------------------


def safe_statistics(
    values,
):
    """
    Compute common statistics.

    Returns
    -------
    tuple
    """

    values = (
        pd.Series(
            values,
            dtype="float64",
        )
        .dropna()
    )

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


def compute_rdf(
    structure: Structure,
    rmax: float,
    bin_size: float,
):
    """
    Compute the radial distribution function using pymatviz.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
    """

    n_bins = max(
        1,
        int(rmax / bin_size),
    )

    r, rdf = calculate_rdf(
        structure=structure,
        cutoff=rmax,
        n_bins=n_bins,
        pbc=(1, 1, 1),
    )

    return (
        np.asarray(r, dtype=float),
        np.asarray(rdf, dtype=float),
    )


# -----------------------------------------------------------------------------


def analyze_structure(
    path: Path,
    rmax: float,
    bin_size: float,
):
    """
    Analyze one crystal structure.

    Parameters
    ----------
    path

    rmax

    bin_size

    Returns
    -------
    tuple
        analysis dictionary
        RDF dataframe
    """

    try:

        structure = load_structure(
            path,
        )

        r, rdf = compute_rdf(
            structure,
            rmax,
            bin_size,
        )

        statistics = extract_rdf_statistics(
            r,
            rdf,
        )

        rdf_frame = pd.DataFrame(
            {
                "file": path.name,
                "radius": r,
                "rdf": rdf,
            }
        )

        result = {
            "file": path.name,
            "formula": structure.composition.formula,
            "reduced_formula": structure.composition.reduced_formula,
            "num_atoms": len(structure),
            **statistics,
            "valid": True,
        }

        return (
            result,
            rdf_frame,
        )

    except Exception as exc:

        result = {
            "file": path.name,
            "formula": None,
            "reduced_formula": None,
            "num_atoms": np.nan,
            "first_peak_position": np.nan,
            "first_peak_height": np.nan,
            "coordination_shell_radius": np.nan,
            "rdf_integral": np.nan,
            "rdf_maximum": np.nan,
            "valid": False,
            "error": str(exc),
        }

        return (
            result,
            None,
        )
# -----------------------------------------------------------------------------


def analyze_dataset(
    files,
    workers,
    rmax,
    bin_size,
):
    """
    Analyze all structures.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
    """

    iterator = tqdm(
        files,
        desc="Computing RDF",
    )

    results = Parallel(
        n_jobs=workers,
    )(
        delayed(
            analyze_structure,
        )(
            path,
            rmax,
            bin_size,
        )
        for path in iterator
    )

    analysis_rows = []
    rdf_frames = []

    for row, rdf in results:

        analysis_rows.append(
            row,
        )

        if rdf is not None:
            rdf_frames.append(
                rdf,
            )

    analysis = pd.DataFrame(
        analysis_rows,
    )

    if rdf_frames:

        rdf_curves = pd.concat(
            rdf_frames,
            ignore_index=True,
        )

    else:

        rdf_curves = pd.DataFrame(
            columns=[
                "file",
                "radius",
                "rdf",
            ]
        )

    return (
        analysis,
        rdf_curves,
    )


# -----------------------------------------------------------------------------


def print_summary(
    df: pd.DataFrame,
):
    """
    Print dataset summary.
    """

    valid = int(
        df["valid"].sum()
    )

    invalid = len(df) - valid

    print("=" * 60)
    print("RDF Analysis Summary")
    print("=" * 60)
    print()

    print(f"Structures            : {len(df)}")
    print(f"Successfully analyzed : {valid}")
    print(f"Failed                : {invalid}")
    print()

    if valid == 0:
        return

    metrics = [
        "first_peak_position",
        "first_peak_height",
        "rdf_integral",
        "rdf_maximum",
    ]

    for column in metrics:

        values = (
            df.loc[
                df["valid"],
                column,
            ]
            .dropna()
        )

        if values.empty:
            continue

        print(column)

        print(
            f"  Mean   : {values.mean():.4f}"
        )

        print(
            f"  Median : {values.median():.4f}"
        )

        print(
            f"  Std    : {values.std():.4f}"
        )

        print(
            f"  Min    : {values.min():.4f}"
        )

        print(
            f"  Max    : {values.max():.4f}"
        )

        print()

    print("Most common formulas")

    counts = Counter(
        df.loc[
            df["valid"],
            "reduced_formula",
        ]
    )

    for formula, count in counts.most_common(10):

        print(
            f"{formula:20s} {count}"
        )

    print()

def extract_rdf_statistics(
    r,
    rdf,
):
    """
    Extract numerical RDF descriptors.
    """

    if len(rdf) == 0:
        return {
            "first_peak_position": np.nan,
            "first_peak_height": np.nan,
            "coordination_shell_radius": np.nan,
            "rdf_integral": np.nan,
            "rdf_maximum": np.nan,
        }

    peaks, _ = find_peaks(
        rdf,
        prominence=0.05,
    )

    if len(peaks):
        peak = peaks[0]

        peak_position = float(r[peak])
        peak_height = float(rdf[peak])

    else:
        peak_position = np.nan
        peak_height = np.nan

    return {
        "first_peak_position": peak_position,
        "first_peak_height": peak_height,
        "coordination_shell_radius": peak_position,
        "rdf_integral": float(trapezoid(rdf, r)),
        "rdf_maximum": float(np.max(rdf)),
    }

# -----------------------------------------------------------------------------


def save_csv(
    df: pd.DataFrame,
    rdf_curves: pd.DataFrame,
    output_dir: Path,
):
    """
    Save CSV outputs.
    """

    df.to_csv(
        output_dir / "rdf_analysis.csv",
        index=False,
    )

    rdf_curves.to_csv(
        output_dir / "rdf_curves.csv",
        index=False,
    )

    summary = [
        {
            "Metric": "Structures",
            "Value": len(df),
        },
        {
            "Metric": "Valid structures",
            "Value": int(
                df["valid"].sum()
            ),
        },
        {
            "Metric": "Invalid structures",
            "Value": int(
                (~df["valid"]).sum()
            ),
        },
    ]

    valid = df[
        df["valid"]
    ]

    numeric_columns = [
        "first_peak_position",
        "first_peak_height",
        "coordination_shell_radius",
        "rdf_integral",
        "rdf_maximum",
    ]

    for column in numeric_columns:

        values = (
            valid[column]
            .dropna()
        )

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


def plot_average_rdf(
    rdf_curves: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot the average RDF over all structures.
    """

    if rdf_curves.empty:
        return

    average = (
        rdf_curves
        .groupby("radius")["rdf"]
        .mean()
        .reset_index()
    )

    fig, ax = plt.subplots()

    ax.plot(
        average["radius"],
        average["rdf"],
        linewidth=2,
    )

    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("g(r)")
    ax.set_title("Average Radial Distribution Function")

    save_figure(
        fig,
        output_dir,
        "average_rdf",
    )


# -----------------------------------------------------------------------------


def plot_rdf_heatmap(
    rdf_curves: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot RDF heatmap.
    """

    if rdf_curves.empty:
        return

    pivot = rdf_curves.pivot(
        index="file",
        columns="radius",
        values="rdf",
    )

    fig, ax = plt.subplots(
        figsize=(10, 8),
    )

    sns.heatmap(
        pivot,
        cmap="viridis",
        ax=ax,
        cbar_kws={
            "label": "g(r)"
        },
    )

    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("Structure")
    ax.set_title("Radial Distribution Functions")

    save_figure(
        fig,
        output_dir,
        "rdf_heatmap",
    )


# -----------------------------------------------------------------------------


def plot_histogram(
    values,
    xlabel: str,
    title: str,
    filename: str,
    output_dir: Path,
    bins: int = 40,
):
    """
    Plot histogram.
    """

    values = pd.Series(
        values,
    ).dropna()

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

    ax.set_xlabel(
        xlabel,
    )

    ax.set_ylabel(
        "Count",
    )

    ax.set_title(
        title,
    )

    save_figure(
        fig,
        output_dir,
        filename,
    )


# -----------------------------------------------------------------------------


def create_all_plots(
    df: pd.DataFrame,
    rdf_curves: pd.DataFrame,
    output_dir: Path,
):
    """
    Generate all figures.
    """

    valid = df[
        df["valid"]
    ]

    if valid.empty:
        return

    plot_average_rdf(
        rdf_curves,
        output_dir,
    )

    plot_rdf_heatmap(
        rdf_curves,
        output_dir,
    )

    plot_histogram(
        valid["first_peak_position"],
        xlabel="Distance (Å)",
        title="First RDF Peak Position",
        filename="first_peak_distribution",
        output_dir=output_dir,
    )

    plot_histogram(
        valid["first_peak_height"],
        xlabel="g(r)",
        title="First RDF Peak Height",
        filename="first_peak_height_distribution",
        output_dir=output_dir,
    )

    plot_histogram(
        valid["coordination_shell_radius"],
        xlabel="Distance (Å)",
        title="First Coordination Shell",
        filename="coordination_shell_distribution",
        output_dir=output_dir,
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

    print(
        f"Found {len(files)} structures."
    )

    if len(files) == 0:

        print(
            "No supported files found."
        )

        return

    print()

    df, rdf_curves = analyze_dataset(
        files=files,
        workers=args.workers,
        rmax=args.rmax,
        bin_size=args.bin_size,
    )

    print_summary(
        df,
    )

    save_csv(
        df,
        rdf_curves,
        output_dir,
    )

    print(
        "Generating plots..."
    )

    create_all_plots(
        df,
        rdf_curves,
        output_dir,
    )

    print()

    print(
        f"Results written to {output_dir.resolve()}"
    )

    print(
        "Done."
    )


# -----------------------------------------------------------------------------


if __name__ == "__main__":
    main()
