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


