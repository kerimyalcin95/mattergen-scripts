#!/usr/bin/env python3
"""
06_fingerprint_analysis.py

Analyze crystal fingerprints for generated crystal structures.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes multiple crystal fingerprints
    * CrystalNNFingerprint
    * OPSiteFingerprint
    * VoronoiFingerprint
- Computes:
    * Mean fingerprint vector
    * Fingerprint norm
    * Pairwise cosine similarity
    * Nearest-neighbor similarity
    * PCA embedding
- Saves:
    * CSV with one row per structure
    * Individual fingerprint matrices
    * Similarity matrix
    * PCA coordinates
    * Summary CSV
    * Failed files CSV
    * Publication-quality PDF and PNG figures
- Prints dataset summary

Supported libraries
-------------------
- pymatgen
- matminer
- ase
- numpy
- pandas
- scipy
- scikit-learn
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

from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

from sklearn.metrics.pairwise import cosine_similarity
from sklearn.decomposition import PCA

from ase.io import read as ase_read

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from matminer.featurizers.site import (
    CrystalNNFingerprint,
    OPSiteFingerprint,
    VoronoiFingerprint,
)

from matminer.featurizers.structure import SiteStatsFingerprint

import warnings

warnings.filterwarnings(
    "ignore",
    message="CrystalNN: cannot locate an appropriate radius.*",
)

# -----------------------------------------------------------------------------

SUPPORTED_EXTENSIONS = {
    ".cif",
    ".xyz",
    ".extxyz",
}

# -----------------------------------------------------------------------------

sns.set_theme(
    style="whitegrid",
)

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
        description="Analyze crystal fingerprints."
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
    Find supported crystal structure files.
    """

    iterator = (
        folder.rglob("*")
        if recursive
        else folder.glob("*")
    )

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

    Returns
    -------
    pymatgen.core.Structure
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
    Create the fingerprint output directory.
    """

    output_dir = (
        base_output
        / "fingerprint"
    )

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
        mean, median, std, min, max
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


def create_featurizers():
    """
    Create structure fingerprint featurizers.
    """

    crystalnn = SiteStatsFingerprint(
        CrystalNNFingerprint.from_preset("ops")
    )

    opsite = SiteStatsFingerprint(
        OPSiteFingerprint()
    )

    voronoi = SiteStatsFingerprint(
        VoronoiFingerprint()
    )

    return (
        crystalnn,
        opsite,
        voronoi,
    )


# -----------------------------------------------------------------------------


def featurize_structure(
    structure: Structure,
    featurizer,
):
    """
    Compute a structure fingerprint.

    Returns
    -------
    numpy.ndarray
    """

    return np.asarray(
        featurizer.featurize(
            structure,
        ),
        dtype=float,
    )


# -----------------------------------------------------------------------------


def compute_fingerprint_norm(
    fingerprint: np.ndarray,
):
    """
    Compute the Euclidean norm of a fingerprint.

    Parameters
    ----------
    fingerprint

    Returns
    -------
    float
    """

    return float(
        np.linalg.norm(
            fingerprint,
        )
    )


# -----------------------------------------------------------------------------


def compute_all_fingerprints(
    structure: Structure,
):
    """
    Compute all fingerprints for a structure.

    Returns
    -------
    dict
    """

    (
        crystalnn,
        opsite,
        voronoi,
    ) = create_featurizers()

    crystalnn_fp = featurize_structure(
        structure,
        crystalnn,
    )

    opsite_fp = featurize_structure(
        structure,
        opsite,
    )

    voronoi_fp = featurize_structure(
        structure,
        voronoi,
    )

    return {
        "crystalnn": crystalnn_fp,
        "opsite": opsite_fp,
        "voronoi": voronoi_fp,
    }


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
    tuple
        analysis dictionary
        fingerprint dictionary
    """

    try:

        structure = load_structure(
            path,
        )

        fingerprints = compute_all_fingerprints(
            structure,
        )

        result = {
            "file": path.name,
            "formula": structure.composition.formula,
            "reduced_formula": structure.composition.reduced_formula,
            "num_atoms": len(structure),
            "crystalnn_norm": compute_fingerprint_norm(
                fingerprints["crystalnn"],
            ),
            "opsite_norm": compute_fingerprint_norm(
                fingerprints["opsite"],
            ),
            "voronoi_norm": compute_fingerprint_norm(
                fingerprints["voronoi"],
            ),
            "valid": True,
        }

        return (
            result,
            fingerprints,
        )

    except Exception as exc:

        return (
            {
                "file": path.name,
                "formula": None,
                "reduced_formula": None,
                "num_atoms": np.nan,
                "crystalnn_norm": np.nan,
                "opsite_norm": np.nan,
                "voronoi_norm": np.nan,
                "valid": False,
                "error": str(exc),
            },
            None,
        )

# -----------------------------------------------------------------------------


def build_fingerprint_dataframe(
    fingerprints,
    key: str,
):
    """
    Convert one fingerprint type into a DataFrame.

    Parameters
    ----------
    fingerprints

    key

    Returns
    -------
    pandas.DataFrame
    """

    rows = []

    for file_name, fp_dict in fingerprints.items():

        vector = fp_dict[key]

        row = {
            "file": file_name,
        }

        row.update(
            {
                f"feature_{i:03d}": value
                for i, value in enumerate(vector)
            }
        )

        rows.append(
            row,
        )

    return pd.DataFrame(
        rows,
    )


# -----------------------------------------------------------------------------


def compute_similarity_statistics(
    matrix: pd.DataFrame,
):
    """
    Compute pairwise cosine similarities.

    Parameters
    ----------
    matrix

    Returns
    -------
    tuple
        similarity dataframe
        nearest-neighbor similarity
        mean similarity
    """

    feature_matrix = (
        matrix
        .drop(columns="file")
        .to_numpy(dtype=float)
    )

    similarity = cosine_similarity(
        feature_matrix,
    )

    similarity_df = pd.DataFrame(
        similarity,
        index=matrix["file"],
        columns=matrix["file"],
    )

    nearest_neighbor = []

    mean_similarity = []

    for row in similarity:

        values = row.copy()

        values = values[
            values < 0.999999
        ]

        if len(values):

            nearest_neighbor.append(
                float(
                    np.max(values)
                )
            )

            mean_similarity.append(
                float(
                    np.mean(values)
                )
            )

        else:

            nearest_neighbor.append(
                np.nan,
            )

            mean_similarity.append(
                np.nan,
            )

    return (
        similarity_df,
        nearest_neighbor,
        mean_similarity,
    )


# -----------------------------------------------------------------------------


def compute_pca(
    matrix: pd.DataFrame,
):
    """
    Compute PCA coordinates.

    Parameters
    ----------
    matrix

    Returns
    -------
    pandas.DataFrame
    """

    feature_matrix = (
        matrix
        .drop(columns="file")
        .to_numpy(dtype=float)
    )

    pca = PCA(
        n_components=2,
        random_state=42,
    )

    coordinates = pca.fit_transform(
        feature_matrix,
    )

    return pd.DataFrame(
        {
            "file": matrix["file"],
            "pc1": coordinates[:, 0],
            "pc2": coordinates[:, 1],
        }
    )


# -----------------------------------------------------------------------------


def analyze_dataset(
    files,
    workers,
):
    """
    Analyze all structures.

    Parameters
    ----------
    files

    workers

    Returns
    -------
    tuple
    """

    iterator = tqdm(
        files,
        desc="Computing fingerprints",
    )

    results = Parallel(
        n_jobs=workers,
    )(
        delayed(
            analyze_structure,
        )(path)
        for path in iterator
    )

    analysis_rows = []

    fingerprint_dict = {}

    for row, fingerprints in results:

        analysis_rows.append(
            row,
        )

        if fingerprints is None:
            continue

        fingerprint_dict[
            row["file"]
        ] = fingerprints

    analysis = pd.DataFrame(
        analysis_rows,
    )

    valid = analysis[
        analysis["valid"]
    ]

    if valid.empty:

        return (
            analysis,
            None,
            None,
            None,
            None,
            None,
        )

    crystalnn_df = build_fingerprint_dataframe(
        fingerprint_dict,
        "crystalnn",
    )

    opsite_df = build_fingerprint_dataframe(
        fingerprint_dict,
        "opsite",
    )

    voronoi_df = build_fingerprint_dataframe(
        fingerprint_dict,
        "voronoi",
    )

    (
        similarity_matrix,
        nearest_similarity,
        mean_similarity,
    ) = compute_similarity_statistics(
        crystalnn_df,
    )

    analysis.loc[
        analysis["valid"],
        "nearest_neighbor_similarity",
    ] = nearest_similarity

    analysis.loc[
        analysis["valid"],
        "mean_similarity",
    ] = mean_similarity

    pca_coordinates = compute_pca(
        crystalnn_df,
    )

    pca_coordinates = pca_coordinates.merge(
        analysis[
            [
                "file",
                "formula",
                "reduced_formula",
            ]
        ],
        on="file",
        how="left",
    )

    return (
        analysis,
        crystalnn_df,
        opsite_df,
        voronoi_df,
        similarity_matrix,
        pca_coordinates,
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
    print("Fingerprint Analysis Summary")
    print("=" * 60)
    print()

    print(f"Structures            : {len(df)}")
    print(f"Successfully analyzed : {valid}")
    print(f"Failed                : {invalid}")
    print()

    if valid == 0:
        return

    metrics = [
        "crystalnn_norm",
        "opsite_norm",
        "voronoi_norm",
        "mean_similarity",
        "nearest_neighbor_similarity",
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


# -----------------------------------------------------------------------------


def save_csv(
    analysis: pd.DataFrame,
    crystalnn: pd.DataFrame,
    opsite: pd.DataFrame,
    voronoi: pd.DataFrame,
    similarity: pd.DataFrame,
    pca: pd.DataFrame,
    output_dir: Path,
):
    """
    Save all CSV outputs.
    """

    analysis.to_csv(
        output_dir / "fingerprint_analysis.csv",
        index=False,
    )

    crystalnn.to_csv(
        output_dir / "crystalnn_fingerprints.csv",
        index=False,
    )

    opsite.to_csv(
        output_dir / "opsite_fingerprints.csv",
        index=False,
    )

    voronoi.to_csv(
        output_dir / "voronoi_fingerprints.csv",
        index=False,
    )

    similarity.to_csv(
        output_dir / "similarity_matrix.csv",
    )

    pca.to_csv(
        output_dir / "pca_coordinates.csv",
        index=False,
    )

    summary = [
        {
            "Metric": "Structures",
            "Value": len(analysis),
        },
        {
            "Metric": "Valid structures",
            "Value": int(
                analysis["valid"].sum()
            ),
        },
        {
            "Metric": "Invalid structures",
            "Value": int(
                (~analysis["valid"]).sum()
            ),
        },
    ]

    valid = analysis[
        analysis["valid"]
    ]

    numeric_columns = [
        "crystalnn_norm",
        "opsite_norm",
        "voronoi_norm",
        "mean_similarity",
        "nearest_neighbor_similarity",
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

    failed = analysis[
        analysis["valid"] == False
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


def plot_similarity_histogram(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot nearest-neighbor similarity histogram.
    """

    values = (
        df["nearest_neighbor_similarity"]
        .dropna()
    )

    if values.empty:
        return

    fig, ax = plt.subplots()

    sns.histplot(
        values,
        bins=40,
        kde=True,
        stat="count",
        ax=ax,
    )

    ax.set_xlabel(
        "Nearest Neighbor Cosine Similarity",
    )

    ax.set_ylabel(
        "Count",
    )

    ax.set_title(
        "Fingerprint Similarity Distribution",
    )

    save_figure(
        fig,
        output_dir,
        "fingerprint_similarity_histogram",
    )


# -----------------------------------------------------------------------------


def plot_similarity_heatmap(
    similarity: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot fingerprint similarity heatmap.
    """

    if similarity.empty:
        return

    fig, ax = plt.subplots(
        figsize=(10, 8),
    )

    sns.heatmap(
        similarity,
        cmap="viridis",
        ax=ax,
        cbar_kws={
            "label": "Cosine Similarity",
        },
    )

    ax.set_xlabel(
        "Structure",
    )

    ax.set_ylabel(
        "Structure",
    )

    ax.set_title(
        "Fingerprint Similarity Matrix",
    )

    save_figure(
        fig,
        output_dir,
        "fingerprint_heatmap",
    )


# -----------------------------------------------------------------------------


def plot_pca(
    pca: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot PCA embedding.
    """

    if pca.empty:
        return

    fig, ax = plt.subplots()

    ax.scatter(
        pca["pc1"],
        pca["pc2"],
        s=30,
    )

    ax.set_xlabel(
        "Principal Component 1",
    )

    ax.set_ylabel(
        "Principal Component 2",
    )

    ax.set_title(
        "Crystal Fingerprint PCA",
    )

    save_figure(
        fig,
        output_dir,
        "fingerprint_pca",
    )


# -----------------------------------------------------------------------------


def plot_norm_boxplot(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot fingerprint norm distributions.
    """

    columns = [
        "crystalnn_norm",
        "opsite_norm",
        "voronoi_norm",
    ]

    plot_df = df[
        columns
    ].copy()

    if plot_df.dropna(
        how="all",
    ).empty:
        return

    fig, ax = plt.subplots()

    sns.boxplot(
        data=plot_df,
        ax=ax,
    )

    ax.set_xticklabels(
        [
            "CrystalNN",
            "OPSite",
            "Voronoi",
        ],
    )

    ax.set_ylabel(
        "Fingerprint Norm",
    )

    ax.set_title(
        "Fingerprint Norm Distribution",
    )

    save_figure(
        fig,
        output_dir,
        "fingerprint_boxplot",
    )


# -----------------------------------------------------------------------------


def create_all_plots(
    analysis: pd.DataFrame,
    similarity: pd.DataFrame,
    pca: pd.DataFrame,
    output_dir: Path,
):
    """
    Generate all figures.
    """

    valid = analysis[
        analysis["valid"]
    ]

    if valid.empty:
        return

    plot_similarity_histogram(
        valid,
        output_dir,
    )

    plot_similarity_heatmap(
        similarity,
        output_dir,
    )

    plot_pca(
        pca,
        output_dir,
    )

    plot_norm_boxplot(
        valid,
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

    print(
        f"Found {len(files)} structures."
    )

    if len(files) == 0:

        print(
            "No supported files found."
        )

        return

    print()

    (
        analysis,
        crystalnn,
        opsite,
        voronoi,
        similarity,
        pca,
    ) = analyze_dataset(
        files,
        args.workers,
    )

    print_summary(
        analysis,
    )

    save_csv(
        analysis,
        crystalnn,
        opsite,
        voronoi,
        similarity,
        pca,
        output_dir,
    )

    print(
        "Generating plots..."
    )

    create_all_plots(
        analysis,
        similarity,
        pca,
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
