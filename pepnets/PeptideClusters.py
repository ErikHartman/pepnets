import pandas as pd

class PeptideClusters:
    def __init__(self, clusters: list()):
        self.clusters = clusters
        self.n_clusters = len(clusters)
        self.proteins = self.get_proteins()

    def reindex(self):
        clusters = self.clusters
        for protein in self.proteins:
            pr_cl = sorted(
                [cluster for cluster in self.clusters if cluster.protein == protein],
                key=lambda x: x.start,
            )
            for i, cluster in enumerate(pr_cl):
                cluster.id = f"{protein}_{i}"
        self.clusters = clusters

    def get_proteins(self):
        proteins = []
        for cluster in self.clusters:
            proteins.append(cluster.protein)
        return sorted(list(set(proteins)))

    def get_cluster_by_id(self, cluster_id: str):
        for cluster in self.clusters:
            if cluster_id == cluster.id:
                return cluster
        print(f"cluster {cluster_id} not found")
        return None

    def add_cluster(
        self,
        cluster,
    ):
        self.clusters.append(cluster)
        self.n_clusters = len(self.clusters)

    def remove_cluster(self, cluster):
        self.clusters.remove(cluster)
        self.n_clusters = len(self.clusters)

    def merge_nearby_clusters(self, wiggle_room=5):
        c = 0
        for cluster1 in self.clusters:
            start1 = cluster1.start
            end1 = cluster1.end
            protein = cluster1.protein
            clusters_in_same_protein = [
                cluster for cluster in self.clusters if protein == cluster.protein
            ]

            for cluster2 in clusters_in_same_protein:
                if cluster1.id is not cluster2.id:
                    if abs(start1 - cluster2.start) <= wiggle_room:
                        if abs(end1 - cluster2.end) <= wiggle_room:
                            for peptide in cluster2.peptides:
                                cluster1.add_peptide(peptide)
                                cluster2.remove_peptide(peptide)
                            self._remove_empty_clusters()
                            c += 1
        print(f"Merged {c} clusters\n")

    def _remove_empty_clusters(self):
        removed_clusters = 0
        new_clusters = []
        for cluster in self.clusters:
            if not cluster.is_empty():
                new_clusters.append(cluster)
            else:
                removed_clusters += 1
        self.clusters = new_clusters
        self.n_clusters = len(self.clusters)
        self.proteins = self.get_proteins()

    def remove_small_clusters(self, threshold):
        for n in range(threshold):
            self.remove_clusters_with_n_peptides(n)
        self.proteins = self.get_proteins()

    def remove_clusters_with_n_peptides(self, n):
        removed_clusters = 0
        new_clusters = []
        for cluster in self.clusters:
            if not cluster.n_peptides == n:
                new_clusters.append(cluster)
            else:
                removed_clusters += 1
        print(f"Removed {removed_clusters} clusters")
        self.clusters = new_clusters
        self.n_clusters = len(self.clusters)
        self.proteins = self.get_proteins()

    def get_n_biggest_clusters(self, n):
        sorted_clusters = sorted(
            self.clusters, key=lambda x: x.n_peptides, reverse=True
        )
        return sorted_clusters[:n]

    def get_cluster(self, peptide, protein):
        for cluster in self.clusters:
            if protein == cluster.protein:
                if peptide in cluster.peptide_sequences:
                    return cluster
        print(f"cluster for {peptide} ({protein}) not found")
        return None

    def to_edgelist(self, savepath):
        with open(savepath, "w") as f:
            f.write(f"from\tto\n")
            for cluster in self.clusters:
                for peptide in cluster.peptide_sequences:
                    f.write(
                        f"{peptide}\t{cluster.id}: ({cluster.start}-{cluster.end})\n"
                    )

                    f.write(
                        f"{cluster.id}: ({cluster.start}-{cluster.end})\t{cluster.protein}\n"
                    )

    def to_feature_edgelist(self, savepath):
        with open(savepath, "w") as f:
            f.write(f"from\tto\n")
            for cluster in self.clusters:
                f.write(
                    f"{cluster.id}_n_peptides\t{cluster.id}: ({cluster.start}-{cluster.end})\n"
                )
                f.write(
                    f"{cluster.id}_intensity\t{cluster.id}: ({cluster.start}-{cluster.end})\n"
                )
                f.write(
                    f"{cluster.id}: ({cluster.start}-{cluster.end})\t{cluster.protein}\n"
                )

    def _get_aa(self, protein, index, database):
        sequence = database[database["Entry Name"] == f"{protein}"]["Sequence"].values[0]
        try:
            return sequence[index]
        except:
            return None
        
     
    def _get_flanks(self, protein, start, end, database):  
        np1 = self._get_aa(protein, start-1, database)
        np1p = self._get_aa(protein, start, database)
        cp1 = self._get_aa(protein, end-1, database)
        cp1p = self._get_aa(protein, end, database)
        return np1, np1p, cp1, cp1p
    
    def to_tsv(self, savepath, database):
        with open(savepath, "w") as f:
            f.write(f"ID\tProtein\tStart\tEnd\tPeptides\tNp1\tNp1p\tCp1\tCp1p\tLongest\n")
            for cluster in self.clusters:
                start = cluster.start
                end = cluster.end
                protein = cluster.protein
                peptide_sequences = cluster.peptide_sequences
                longest_peptide = sorted(peptide_sequences, key=lambda x: len(x), reverse=True)[0]
                np1, np1p, cp1, cp1p = self._get_flanks(protein, start, end, database)
                f.write(
                    f"{cluster.id}\t{protein}\t{start}\t{end}\t{peptide_sequences}\t{np1}\t{np1p}\t{cp1}\t{cp1p}\t{longest_peptide}\n"
                )

    def to_df(self):
        save_dict = {
            "cluster": [],
            "protein": [],
            "peptide": [],
            "start": [],
            "end": [],
            "n_peptides": [],
        }
        for cluster in self.clusters:
            for peptide in cluster.peptides:
                save_dict["cluster"].append(cluster.id)
                save_dict["protein"].append(peptide.protein)
                save_dict["peptide"].append(peptide.sequence)
                save_dict["start"].append(peptide.start)
                save_dict["end"].append(peptide.end)
                save_dict["n_peptides"].append(cluster.n_peptides)

        return pd.DataFrame(save_dict)

    def __repr__(self) -> str:
        return f"{[cluster.id for cluster in self.clusters]}"

    def __iter__(self):
        return self.clusters.__iter__()
