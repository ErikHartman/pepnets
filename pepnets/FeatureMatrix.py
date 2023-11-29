from pepnets.PeptideClusters import PeptideClusters
import pandas as pd
import numpy as np


class FeatureMatrix:
    def __init__(
        self,
        datamatrix: pd.DataFrame,
        clusters: PeptideClusters,
    ):
        self.datamatrix = datamatrix.replace(0, np.nan)
        self.clusters = clusters
        self.samples = [col for col in datamatrix.columns if "Sample" in col]
        self.peptides = datamatrix["Peptide"].values.tolist()
        self.proteins = datamatrix["Protein"].values.tolist()
        cluster_column = []
        for protein, peptide in zip(self.proteins, self.peptides):
            cluster = clusters.get_cluster(peptide, protein)
            if cluster:
                cluster_column.append(cluster.id)
            else:
                cluster_column.append(np.nan)

        self.datamatrix["Cluster"] = cluster_column
        pre_drop = len(self.datamatrix.index)
        self.datamatrix = self.datamatrix.dropna(subset=["Cluster"])
        print("Dropped: ", pre_drop - len(self.datamatrix.index))

    def get_topn(self):
        datamatrix = self.datamatrix.copy()
        datamatrix = datamatrix.drop(
            columns=["Protein", "Start", "End", "Peptide", "Mods", "id"]
        )
        clusters = self.datamatrix["Cluster"].unique()

        quantifications = {}

        for cluster in clusters:
            dm = (
                datamatrix[datamatrix["Cluster"] == cluster]
                .copy()
                .set_index(["Cluster"])
            )
            quantification = quantify_group(dm.values)

            quantifications[cluster] = quantification

        datamatrix = datamatrix.set_index(["Cluster"])
        q = pd.DataFrame(quantifications).T
        q.columns = datamatrix.columns
        return q


def quantify_group(
    group_data: np.ndarray, summarization_method: str = "sum"
) -> np.ndarray:
    top_n = 3
    group_data = np.nan_to_num(group_data, nan=0.0)

    sort_indices = np.argsort(group_data, axis=0)[::-1]

    sorted_precursors = np.take_along_axis(group_data, sort_indices, axis=0)

    if summarization_method == "sum":
        quantification: np.ndarray = np.sum(sorted_precursors[:top_n], axis=0)

    elif summarization_method == "mean":
        quantification: np.ndarray = np.mean(sorted_precursors[:top_n], axis=0)

    elif summarization_method == "median":
        quantification: np.ndarray = np.median(sorted_precursors[:top_n], axis=0)

    return quantification
