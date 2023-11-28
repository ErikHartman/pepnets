import numpy as np
from collections import Counter


class PeptideCluster:
    def __init__(self, cluster_id, peptides, protein):
        self.id = cluster_id
        self.peptides = peptides
        self.protein = protein
        self.n_peptides = len(peptides)
        self.peptide_sequences = self._get_unique_peptides()
        self.start, self.end = self._get_endpoints()
        self.inter_cluster_distance = self._get_inter_cluster_distance()

    def _get_unique_peptides(self):
        peptides = set()
        for peptide in self.peptides:
            peptides.add(peptide.sequence)
        return list(peptides)

    def _get_endpoints(self, method="mode"):
        starts = []
        ends = []
        for peptide in self.peptides:
            starts.append(peptide.start)
            ends.append(peptide.end)
        if len(self.peptides) == 0:
            return 0, 0
        if method == "mode":
            start_counter = Counter(starts)
            end_counter = Counter(ends)
            if max(start_counter.values()) == 1:
                start = min(starts)
            else:
                start = max(set(starts), key=starts.count)
            if max(end_counter.values()) == 1:
                end = max(ends)
            else:
                end = max(set(ends), key=ends.count)
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

    def __repr__(self) -> str:
        return f"{self.id}, proteins: ({self.protein}), #peptides: {self.n_peptides}"
