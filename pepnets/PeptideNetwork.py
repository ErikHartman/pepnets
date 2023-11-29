import igraph as ig
import pandas as pd
import networkx as nx
import leidenalg
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm

from pepnets.PeptideCluster import PeptideCluster
from pepnets.PeptideClusters import PeptideClusters
from pepnets.Peptide import Peptide



plt.rcParams["font.family"] = "Arial"


class PeptideNetwork:
    def __init__(
        self,
        datamatrix: pd.DataFrame,
        protein_database: pd.DataFrame,
    ):
        self.datamatrix = datamatrix.sort_values("Start", ascending=False)
        self.protein_database = protein_database
        if protein_database is not None:
            proteins_in_database = protein_database["Entry Name"].values.tolist()
            self.datamatrix = self.datamatrix[
                self.datamatrix["Protein"].isin(proteins_in_database)
            ]
        self.peptides, self.proteins = self._generate_peptides()

    def _generate_peptides(self):
        print("Reading peptides...")
        peptides = []
        proteins = []
        id = 0
        for _, row in self.datamatrix.iterrows():
            protein = row["Protein"]
            peptide = row["Peptide"]
            start = self._get_peptide_start(peptide, row["Protein"])
            peptides.append(Peptide(peptide, start, protein, id))
            proteins.append(protein)
            id += 1
        proteins = list(set(proteins))
        return peptides, proteins


    def create_network(self, distance_cutoff: int = 4):
        protein_networks = {}
        for protein in self.proteins:
            peptides_in_protein = [
                peptide for peptide in self.peptides if peptide.protein == protein
            ]
            G = nx.Graph()
            for peptide in peptides_in_protein:
                G.add_node(
                    peptide.id,
                    peptide=peptide.sequence,
                    protein=peptide.protein,
                    start=peptide.start,
                    end=peptide.end,
                ),

            for i, peptide_from in enumerate(peptides_in_protein):
                for peptide_to_index in range(i + 1):
                    peptide_to = peptides_in_protein[peptide_to_index]

                    overlap_percentage = get_overlap_percentage(
                        peptide_from, peptide_to
                    )
                    epsilon = 1e-8

                    if peptide_from.sequence == peptide_to.sequence:
                        continue
                    if overlap_percentage == 0:
                        continue

                    overlap_distance = 1 / (overlap_percentage + epsilon)

                    length_ratio = get_length_ratio(peptide_from, peptide_to)

                    center_distance = get_center_distance(peptide_from, peptide_to)

                    d = length_ratio + overlap_distance + center_distance

                    if d <= distance_cutoff:
                        d = d + epsilon
                        G.add_edge(
                            peptide_from.id,
                            peptide_to.id,
                            distance=d,
                            distance_inv=1 / d,
                        )
            protein_networks[protein] = G
        self.protein_networks = protein_networks

        return protein_networks

    def get_clusters(self, resolution: float, random_seed: int = 42):
        peptide_clusters = []
        for protein, G in zip(
            self.protein_networks.keys(), self.protein_networks.values()
        ):
            H = ig.Graph.from_networkx(G)
            partition = leidenalg.find_partition(
                H,
                leidenalg.RBConfigurationVertexPartition,
                n_iterations=4,
                weights="distance_inv",
                seed=random_seed,
                resolution_parameter=resolution,
            )
            clusters = {}
            for i, cluster in enumerate(partition):
                for node in cluster:
                    peptide_id = H.vs[node]["_nx_name"]
                    clusters[peptide_id] = i

            cluster_dict = {}
            for cluster_index, peptide_index in zip(clusters.values(), clusters.keys()):
                cluster_name = f"{protein}_{cluster_index}"
                if cluster_name not in cluster_dict.keys():
                    cluster_dict[cluster_name] = []
                peptide = self._get_peptide(peptide_index)
                cluster_dict[cluster_name].append(peptide)

            for cluster_id in cluster_dict.keys():
                peptide_clusters.append(
                    PeptideCluster(
                        cluster_id=cluster_id,
                        peptides=cluster_dict[cluster_id],
                        protein=protein,
                    )
                )
        clusters = PeptideClusters(peptide_clusters)
        self.clusters = clusters
        return clusters

    def _get_peptide(self, id):
        for peptide in self.peptides:
            if peptide.id == id:
                return peptide
        return None

    def _get_protein_sequence(self, protein):
        if protein in self.protein_database["Entry Name"].values.tolist():
            protein_seq = self.protein_database[
                self.protein_database["Entry Name"] == protein
            ]["Sequence"].values[0]
            return protein_seq
        return "Not in database"

    def _get_peptide_start(self, peptide, protein):
        if self.protein_database is not None:
            if protein in self.protein_database["Entry Name"].values.tolist():
                protein_seq = self.protein_database[
                    self.protein_database["Entry Name"] == protein
                ]["Sequence"].values[0]
                start = protein_seq.find(peptide)
                return start
            print(f"Protein {protein} not in database.")
        else:
            return self.datamatrix[self.datamatrix["Peptide"] == peptide][
                "Start"
            ].values[0]
        return "Not in database"

    def plot_protein(self, protein, save_str=None, figsize=(8,4)):
        plt.clf()
        fig, axs = plt.subplots(1, 2, figsize=figsize)
        graph = self.protein_networks[protein]
        H = ig.Graph.from_networkx(graph)
        partition = leidenalg.find_partition(
            H,
            leidenalg.RBConfigurationVertexPartition,
            n_iterations=4,
            weights="distance_inv",
            resolution_parameter=.8,
        )
        clusters = {}
        for i, cluster in enumerate(partition):
            for node in cluster:
                peptide_id = H.vs[node]["_nx_name"]
                clusters[peptide_id] = i
        partition=clusters
        layout = nx.spring_layout(
            graph, weight="distance_inv", k=1.5 / np.sqrt(graph.number_of_nodes())
        )

        cmap_partition = cm.get_cmap("tab20", max(partition.values()) + 1)

        n1 = nx.draw_networkx_nodes(
            graph,
            layout,
            partition.keys(),
            node_size=1,
            node_color=[graph.nodes[n]["start"] for n in graph.nodes()],
            cmap=plt.cm.Blues,
            alpha=0.8,
            ax=axs[0],
        )
        plt.colorbar(n1, ax=axs[0], label="Start position")
        nx.draw_networkx_edges(graph, layout, alpha=0.5, ax=axs[0])

        n2 = nx.draw_networkx_nodes(
            graph,
            layout,
            partition.keys(),
            node_size=3,
            node_color=list(partition.values()),
            cmap=cmap_partition,
            alpha=0.8,
            ax=axs[1],
        )
        plt.colorbar(n2, ax=axs[1], label="Designated cluster")
        nx.draw_networkx_edges(graph, layout, alpha=0.5, ax=axs[1])
        plt.suptitle(f"{protein} network")
        sns.despine()
        plt.tight_layout()
        if save_str:
            plt.savefig(save_str, dpi=1200)


def get_overlap_percentage(peptide1, peptide2, divisor="total_length"):
    if peptide1.start > peptide2.start:
        peptide1, peptide2 = peptide2, peptide1

    if peptide1.end < peptide2.start:
        return 0

    if peptide2.start >= peptide1.start and peptide2.end <= peptide1.end:
        return 1

    overlap = peptide1.end - peptide2.start
    total_length = peptide1.length + peptide2.length - overlap

    if divisor == "total_length":
        return overlap / (total_length)
    elif divisor == "longest_peptide":
        return overlap / max(peptide1.length, peptide2.length)
    elif divisor == "shortest_peptide":
        return overlap / min(peptide1.length, peptide2.length)
    elif divisor == None:
        return total_length - overlap


def get_center_distance(peptide1, peptide2):
    center_distance = np.abs(peptide1.center - peptide2.center)
    return center_distance


def get_length_difference(peptide1, peptide2):
    return np.abs(peptide1.length - peptide2.length)


def get_length_ratio(peptide1, peptide2):
    return max(peptide1.length / peptide2.length, peptide2.length / peptide1.length)


def get_endpoints_distance(peptide1, peptide2):
    start_dist = np.abs(peptide1.start - peptide2.start)
    end_dist = np.abs(peptide1.end - peptide2.end)
    return start_dist + end_dist
