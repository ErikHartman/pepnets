import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pepnets.palette import *
import seaborn as sns
import logomaker
from pepnets.util import custom_mode


plt.rcParams["font.family"] = "Arial"

amino_acids = list("ARNDBCEQZGHILKMFPSTWYV*")

sns.color_palette("BuPu")

positions = ["p4", "p3", "p2", "p1", "p1'", "p2'", "p3'", "p4'"]
color_scheme = {aa: "#e3e3e3" for aa in amino_acids}
color_scheme.update(
    {
        "K": group_palette["P. aeruginosa"],
        "D": group_palette["S. aureus"],
        "E": group_palette["S. aureus"],
        "S": group_palette["Double infection"],
        "V": sns.color_palette("BuPu")[3],
        "I": sns.color_palette("BuPu")[3],
        "A": sns.color_palette("BuPu")[3],
        "T": sns.color_palette("BuPu")[3],
    }
)


class LogoPlot:
    def __init__(
        self,
        datamatrix: pd.DataFrame,
        database: dict,  # dict of identifier:sequence
        drop: list = None,
        entity : str = "Cluster"
    ):
        self.dm = datamatrix
        self.database = database
        self.drop = drop
        self.entity = entity

    def get_letter_heights(self, test_samples: list, background_samples: list):
        test_nc = self.get_terminal_amino_acid_frequencies(test_samples)
        back_nc = self.get_terminal_amino_acid_frequencies(background_samples)
        height = get_kl(test_nc.T, back_nc.T)

        height["index"] = [-4, -3, -2, -1, 1, 2, 3, 4]
        height.set_index("index", inplace=True)
        return height

    def plot(self, height, ax):
        _ = logomaker.Logo(height, color_scheme=color_scheme, ax=ax)

    def get_terminal_amino_acid_frequencies(self, samples: list):
        nc = self._get_terminal_amino_acid_counts(samples)
        nc_freq = self._get_frequency_df(nc)
        nc = nc_freq / nc_freq.sum(axis=0)
        return nc

    def _get_terminal_amino_acid_counts(self, samples: list,):
        dm = self.dm.copy()
        dm["sum_int"] = dm[samples].sum(axis=1)
        dm = dm[dm["sum_int"] > 0].copy()
        dm = dm[[self.entity, "Protein", "Start", "End", "sum_int"]]
        dm = dm.groupby([self.entity, "Protein"], as_index=False).agg(
            {
                "Start": lambda x: (custom_mode(x, equal_strategy="min")),
                "End": lambda x: (custom_mode(x, equal_strategy="max")),
                "sum_int": "mean",
            }
        )
        
        dm = dm[[self.entity, "Start", "Protein", "End", "sum_int"]]

        n_term = self._get_flanks(dm, term="n")
        c_term = self._get_flanks(dm, term="c")
        nc = pd.concat([n_term, c_term])
        nc.drop(columns={"Start", "End", self.entity}, inplace=True)
        if self.drop:
            nc = n_term[~n_term["p1"].isin(self.drop)]

        return nc

    def _get_frequency_df(self, nc, impute=1e-4) -> pd.DataFrame:
        frequency_dict = {}

        for position in positions:
            frequency_dict[position] = {aa: impute for aa in amino_acids}
        nc["count"] = 1
        for _, row in nc.iterrows():
            value = row["sum_int"]
            for position in positions:
                aa_at_position = row[position]
                if aa_at_position not in ["X", None]:
                    frequency_dict[position][aa_at_position] += value
        return pd.DataFrame(frequency_dict)

    def _get_aa(self, protein, index):
        sequence = self.database[protein]
        try:
            return sequence[index]
        except:
            return None

    def _get_flanks(self, df, term="n"):
        if term == "n":
            start_or_end = "Start"
        else:
            start_or_end = "End"
        df = df.copy()
        
        df["p4"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) - 4),
            axis=1,
        )
        df["p3"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) - 3),
            axis=1,
        )
        df["p2"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) - 2),
            axis=1,
        )
        df["p1"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) - 1),
            axis=1,
        )

        df["p1'"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end])),
            axis=1,
        )
        df["p2'"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) + 1 ),
            axis=1,
        )
        df["p3'"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) + 2),
            axis=1,
        )
        df["p4'"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) + 3 ),
            axis=1,
        )
        return df


def get_kl(p, q):
    """
    KL-divergence
    """
    I = (p * np.log2(p / q)).sum(axis=1)
    return (p.T * I).T
