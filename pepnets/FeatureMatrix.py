from pepnets.PeptideClusters import PeptideClusters
import networkx as nx
import pandas as pd
import numpy as np
import dpks


class FeatureMatrix:
    def __init__(
        self,
        datamatrix: pd.DataFrame,
        design: pd.DataFrame,
        clusters: PeptideClusters,
        min_nr_samples: int = 5,
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
                print(peptide, protein)
                cluster_column.append(np.nan)

        self.datamatrix["Cluster"] = cluster_column
        pre_drop = len(self.datamatrix.index)
        self.datamatrix = self.datamatrix.dropna(subset=["Cluster"])

        self.datamatrix = self.datamatrix.dropna(
            subset=self.samples, thresh=min_nr_samples
        ).fillna(0)
        print("Dropped: ", pre_drop - len(self.datamatrix.index))
        self.design = design
        self.qm = dpks.QuantMatrix(self.datamatrix, self.design)

    def generate_feature_matrix(self):
        maxlfq = (
            self.qm.quantify("maxlfq", level="Cluster").to_df().set_index("Cluster")
        )
        mean = self.datamatrix[self.samples + ["Cluster"]].groupby("Cluster").mean()
        std = self.datamatrix[self.samples + ["Cluster"]].groupby("Cluster").std()
        topn = self.qm.quantify(method="top_n", top_n=3, level="Cluster").to_df().set_index("Cluster")
        protein = self.qm.quantify("top_n", level="Protein").to_df().set_index("Protein")
        n_peptides = (
            self.datamatrix[self.samples + ["Cluster"]]
            .replace(0, np.nan)
            .groupby("Cluster")
            .count()
        )
        return {
            "protein": protein,
            "maxlfq": maxlfq,
            "mean": mean,
            "std": std,
            "n_pep": n_peptides,
            "topn": topn,
        }

