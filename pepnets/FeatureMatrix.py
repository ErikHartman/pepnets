from pepnets.PeptideClusters import PeptideClusters
from pepnets.PeptideCluster import PeptideCluster
import pandas as pd
import numpy as np
from dpks import QuantMatrix
from scipy import stats


class FeatureMatrix:
    def __init__(
        self,
        datamatrix: pd.DataFrame,
        design: pd.DataFrame,
        clusters: PeptideClusters,
    ):
        self.datamatrix = datamatrix
        self.clusters = clusters
        self.design = design
        self.samples = design["sample"].values.tolist()

        cluster_column = []
        cluster_start_column = []
        cluster_end_column = []
        for protein, peptide in zip(
            datamatrix["Protein"].values.tolist(), datamatrix["Peptide"].values.tolist()
        ):
            cluster: PeptideCluster = clusters.get_cluster(peptide, protein)
            if cluster:
                cluster_column.append(cluster.id)
                cluster_start_column.append(cluster.start)
                cluster_end_column.append(cluster.end)
            else:
                cluster_column.append(np.nan)
                cluster_start_column.append(np.nan)
                cluster_end_column.append(np.nan)

        self.datamatrix["Cluster"] = cluster_column
        self.datamatrix["Cluster_start"] = cluster_start_column
        self.datamatrix["Cluster_end"] = cluster_end_column
        self.datamatrix = self.datamatrix.dropna(subset=["Cluster"])

    def get_normalized_dm(self):
        qm = QuantMatrix(self.datamatrix, self.design)
        qm.normalize(
            method="mean",
            use_rt_sliding_window_filter=False,
            log_transform=False,
        )
        return qm.to_df()

    def get_topn(self, normalize=True, scale=False, topn=3):
        qm = QuantMatrix(self.datamatrix, self.design)
        if normalize:
            qm.normalize(
                method="mean",
                use_rt_sliding_window_filter=False,
                log_transform=False,
            )
        if scale:
            qm.scale(method="zscore")
        qm.quantify(
            method="top_n", top_n=topn, summarization_method="mean", level="Cluster"
        )
        return qm.to_df()

    def get_topn_de(
        self,
        group_a,
        group_b,
        min_samples: int = 3,
        top_n=3,
        metric="de",
        normalize=True,
    ):
        """The top N are evaluated based on their level of differential expression."""
        qm = QuantMatrix(self.datamatrix, self.design)
        if normalize:
            qm.normalize(
                method="mean",
                use_rt_sliding_window_filter=False,
                log_transform=False,
            )

        design = self.design.copy()

        data_test = qm.to_df().set_index("Peptide")

        group_1_samples = design[design["group"] == group_a]["sample"]
        group_2_samples = design[design["group"] == group_b]["sample"]

        log_pvals = []
        log_fcs = []

        for _, row in data_test.iterrows():
            row = row[design["sample"]].copy()
            group_1_data = row[group_1_samples]
            group_2_data = row[group_2_samples]

            group_1_n = group_1_data.count()
            group_2_n = group_2_data.count()

            if (group_1_n <= min_samples) or (group_2_n <= min_samples):

                log_pvals.append(np.nan)
                log_fcs.append(np.nan)
                continue

            group_1_data = group_1_data.values.astype(float)
            group_2_data = group_2_data.values.astype(float)

            group_1_data = group_1_data[~np.isnan(group_1_data)]
            group_2_data = group_2_data[~np.isnan(group_2_data)]

            pvalue = stats.ttest_ind(group_1_data, group_2_data).pvalue
            log_pvalue = -np.log10(pvalue)

            group_1_mean = np.nanmean(group_1_data)
            group_2_mean = np.nanmean(group_2_data)

            log_fc = group_1_mean - group_2_mean

            log_pvals.append(log_pvalue)
            log_fcs.append(log_fc)

        data_test["pval"] = log_pvals
        data_test["fc"] = log_fcs

        de_scores = [
            np.sqrt((p / np.nanmax(log_pvals)) ** 2 + (fc / np.nanmax(log_fcs)) ** 2)
            for p, fc in zip(log_pvals, log_fcs)
        ]

        data_test["de"] = de_scores

        unique_clusters = data_test["Cluster"].unique()

        top_n_de = {}

        for cluster in unique_clusters:
            cluster_data = data_test[data_test["Cluster"] == cluster].sort_values(
                metric, ascending=False
            )

            top_de_peptides = cluster_data[0:top_n]

            top_de_peptides = top_de_peptides[design["sample"]]

            samples = top_de_peptides.columns
            top_de_peptides = np.nanmean(top_de_peptides.values, axis=0)

            top_n_de[cluster] = {s: p for s, p in zip(samples, top_de_peptides)}

        return pd.DataFrame(top_n_de).T.reset_index(names=["Cluster"])
