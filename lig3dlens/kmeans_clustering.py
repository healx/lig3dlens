import numpy as np
from sklearn.cluster import MiniBatchKMeans
from scipy.spatial.distance import cdist
import pandas as pd

from typing import Union, List
import click
from loguru import logger

import datamol as dm
from molfeat.calc import FPCalculator
from molfeat.trans import MoleculeTransformer


def find_cluster_centers_opt(
    df: pd.DataFrame, centers: Union[List[List[float]], np.ndarray]
) -> np.ndarray:
    """
    Find molecular cluster centers in a pandas dataframe

    Parameters
    ----------
    df: pd.DataFrame
        dataframe containing compounds and the clusters they are into.
    centers: np.ndarray
        a list or array of cluster centers coordinates.

    Returns
    -------
    np.ndarray: a numpy array representing whether each row in 'df'
               is a cluster center or not.
               The values will be "Yes" for cluster centers
               and "No" for others.
    """

    center_indices = np.argmin(
        cdist(centers, df.fp.tolist(), metric="euclidean"), axis=0
    )
    center_names = df.values[center_indices, 1]
    return np.where(np.isin(df.values[:, 1], center_names), "Yes", "No")


def kmeans_cluster(df: pd.DataFrame, num_clusters: int) -> pd.DataFrame:
    """
    Use kmeans algorithm to cluster the cmpds in the pandas dataframe
    and find the clusters' centers

    Parameters
    ----------
    df: pd.Dataframe
        dataframe containing library compounds
    num_enters: int
        total number of clusters

    Returns
    -------
    df: pd.DataFrame
        df + (clusters and their centers)
    """

    num_rows, _ = df.shape

    # Check if the number of compounds in the df is larger than 10_000
    if num_rows > 10000:
        # number of samples needs to be at least equal the number of clusters
        rows_to_sample = max(int(num_rows / 10), num_clusters)
        train_df = df.sample(rows_to_sample)
        logger.warning("Input file should have less than 10_000 cmpds!")
        logger.warning(f"Sampled {rows_to_sample} rows to train kmeans")
    else:
        train_df = df

    # Numpy array holding the fingerprints of the compounds in the training set
    fp_arr = np.array(train_df.fp.tolist(), dtype=np.int8)

    km = MiniBatchKMeans(
        n_clusters=num_clusters,
        random_state=42,
        batch_size=3 * num_clusters,
        n_init="auto",
    )

    # Fit kmeans to training data
    km.fit(fp_arr)

    # Numpy array with the fingerprints of all the compounds in the input file
    all_data = np.array(df.fp.tolist(), dtype=np.int8)

    # Store cluster ID to df column 'Cluster'
    df["Cluster"] = km.predict(all_data)

    # Store information about cluster centers to df column 'Center'
    df["Center"] = find_cluster_centers_opt(df, km.cluster_centers_)

    return df


def main(infile_name, clusters, outfile_name, fp_dim, fp_type):
    logger.info(f"Loading cmpds from {infile_name} to a Pandas dataframe")

    vs_hits = dm.read_sdf(infile_name, as_df=True, mol_column="ROMol")

    logger.debug(f"There are {vs_hits.shape[0]} cmpds in {infile_name}")
    logger.info(
        f"Generating the {fp_type}::{fp_dim} fingerprints of the provided cmpds"
    )

    # Initialising the fingerprint calculator
    calc = FPCalculator(fp_type, length=fp_dim)

    # Featurisation pipeline
    mol_transf = MoleculeTransformer(calc, n_jobs=-1)

    # Store the fingerprints to the pandas dataframe
    vs_hits["fp"] = mol_transf(vs_hits.ROMol)

    # logger.info(f'Clustering the compounds')
    df_clustered = kmeans_cluster(vs_hits, clusters)

    dm.to_sdf(df_clustered, outfile_name)
    logger.success(f"Library compounds and their clusters stored in {outfile_name}")


if __name__ == "__main__":

    @click.command("input_output")
    @click.option("--in", "infile_name", required=True, help="Input file name")
    @click.option(
        "--out",
        "outfile_name",
        required=True,
        help="Output csv file with SMILES, molecule name, and cluster id",
    )
    @click.option(
        "--clusters",
        "clusters",
        type=int,
        required=True,
        help="Number of clusters to output",
    )
    @click.option(
        "--dim", "fp_dim", type=int, default=1024, help="Number of fingerprint bits"
    )
    @click.option(
        "--fp_type",
        "fp_type",
        default="ecfp",
        type=click.Choice(
            [
                "maccs",
                "ecfp",
                "fcfp",
                "topological",
                "atompair",
                "rdkit",
                "pattern",
                "layered",
                "erg",
                "estate",
                "avalon-count",
                "rdkit-count",
                "ecfp-count",
                "fcfp-count",
                "topological-count",
                "atompair-count",
            ]
        ),
        help=(
                "Fingerprint type (must be one of maccs, ecfp, fcfp, topological, "
                "atompair, rdkit, pattern, layered, erg, estate, avalon-count, "
                "rdkit-count, ecfp-count, fcfp-count, topological-count, "
                "atompair-count)",
        )
    )
    def cli(infile_name, clusters, outfile_name, fp_dim, fp_type):
        main(infile_name, clusters, outfile_name, fp_dim, fp_type)

    cli()
