
class Peptide:
    def __init__(self, sequence: str, start: int, protein: str, id: int):
        self.id = id
        self.sequence = sequence
        self.start = int(start)
        self.protein = protein

        self.length = len(sequence)
        self.end = start + self.length
        self.center = self.start + self.length / 2

    def __repr__(self) -> str:
        return f"{self.sequence} ({self.protein}), start: {self.start}, end: {self.end}"
