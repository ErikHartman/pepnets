import numpy as np
from pepnets.util import custom_mode


class PeptideCluster:
    def __init__(self, cluster_id, peptides, protein):
        self.id = cluster_id
        self.peptides = peptides
        self.protein = protein
        self.n_peptides = len(peptides)
        self.peptide_sequences = self._get_unique_peptides()
        self.start, self.end = self._get_endpoints(method="mode")
        self.inter_cluster_distance = self._get_inter_cluster_distance()

    def _get_unique_peptides(self):
        peptides = set()
        for peptide in self.peptides:
            peptides.add(peptide.sequence)
        return sorted(list(peptides))

    def _get_endpoints(self, method="mode"):
        starts = []
        ends = []
        for peptide in self.peptides:
            starts.append(peptide.start)
            ends.append(peptide.end)
        if len(self.peptides) == 0:
            return 0, 0
        if method == "mode":
            start = custom_mode(starts, equal_strategy="min")
            end = custom_mode(ends, equal_strategy="max")
            return start, end
        elif method == "longest":
            return min(starts), max(ends)

    def add_peptide(self, peptide):
        self.peptides.append(peptide)
        self.n_peptides = len(self.peptides)
        self.peptide_sequences = self._get_unique_peptides()
        self.start, self.end = self._get_endpoints()

    def remove_peptide(self, peptide_to_remove):
        new_peptides = [
            peptide
            for peptide in self.peptides
            if peptide.sequence != peptide_to_remove.sequence
        ]
        self.peptides = new_peptides
        self.n_peptides = len(self.peptides)
        self.peptide_sequences = self._get_unique_peptides()
        self.start, self.end = self._get_endpoints()

    def _get_inter_cluster_distance(self):
        center_distance = 0
        n_distances = 0
        for peptide1 in self.peptides:
            for peptide2 in self.peptides:
                center_distance += np.abs(peptide1.center - peptide2.center)
                n_distances += 1
        return center_distance / n_distances


    def is_empty(self):
        if len(self.peptides) == 0:
            return True
        return False

    def hash(self):
        hash_peptide = ""
        for peptide in self.peptide_sequences:
            hash_peptide += peptide
        return hash_peptide
    
    def name(self):
        return f"{self.protein} ({self.start}-{self.end})"

    def __repr__(self) -> str:
        return f"{self.id}, proteins: ({self.protein}), #peptides: {self.n_peptides}"
