import numpy as np
import pandas as pd
from pepnets.palette import *
import seaborn as sns
import logomaker
from collections import Counter


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
        design: pd.DataFrame,
        database: pd.DataFrame,
        topn: pd.DataFrame,
        species: str = "PIG",
        weight: str = "mean_int",
        drop: list = None
    ):
        self.dm = datamatrix
        self.design = design
        self.database = database
        self.species = species
        self.topn = topn
        self.weight = weight
        self.drop = drop

    def get_letter_heights(self, test_samples: list, background_samples: list):
        test_nc = self.get_nc(test_samples)
        back_nc = self.get_nc(background_samples)
        height = get_height_from_norm(test_nc, back_nc)
        height["index"] = [-4, -3, -2, -1, 1, 2, 3, 4]
        height.set_index("index", inplace=True)
        return height

    def plot(self, height, ax):
        logo = logomaker.Logo(height, color_scheme=color_scheme, ax=ax)

    def get_nc(self, samples: list):
        n, c = self._prepare_data(samples)
        nc = pd.concat([n, c])
        nc_freq = self.get_frequency_dict(nc)
        if self.drop:
            nc_freq = nc_freq.drop(index=self.drop)
        nc = nc_freq / nc_freq.sum(axis=0)
        return nc

    def _prepare_data(self, samples: list) -> (pd.DataFrame, pd.DataFrame):
        """
        Return two dataframes for the n-term and c-term respectively
        """
        dm = self.dm.copy()
        dm["mean_int"] = dm[samples].mean(axis=1)
        dm = dm[dm["mean_int"] > 0].copy()
        dm = dm[["Cluster", "Protein", "Start", "End", "mean_int"]]
        dm = dm.groupby(["Cluster", "Protein"], as_index=False).agg(
            {
                "Start": lambda x:  min(x) if max(Counter(x).values()) == 1 else max(list(x), key=list(x).count),
                "End": lambda x:  max(x) if max(Counter(x).values()) == 1 else max(list(x), key=list(x).count),
                "mean_int": "mean",
            }
        )
        dm = dm[["Cluster", "Start", "Protein", "End", "mean_int"]]
        n_term = self._get_flanks(dm, term="n")
        c_term = self._get_flanks(dm, term="c")
        n_term.drop(columns={"Start", "End", "Cluster"}, inplace=True)
        c_term.drop(columns={"Start", "End", "Cluster"}, inplace=True)

        return n_term, c_term

    def get_frequency_dict(self, nc, impute=0.0001) -> pd.DataFrame:
        """
        Returns a frequency dict from a frequency dataframe
        """
        frequency_dict = {}

        for position in positions:
            frequency_dict[position] = {aa: impute for aa in amino_acids}

        nc["count"] = 1
        for _, row in nc.iterrows():
            value = row[self.weight]
            for position in positions:
                aa_at_position = row[position]
                if aa_at_position not in ["X", None]:
                    frequency_dict[position][aa_at_position] += value
        return pd.DataFrame(frequency_dict)

    def _get_aa(self, protein, index):
        sequence = self.database[
            self.database["Entry Name"].isin([f"{protein}_PIG", f"{protein}_HUMAN"])
        ]["Sequence"].values[0]
        try:
            return sequence[index]
        except:
            return None

    def _get_flanks(self, df, term="n"):
        if term == "n":
            start_or_end = "Start"
            add = 0
        else:
            start_or_end = "End"
            add = 1
        df = df.copy()
        df["p4"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) - 5 + add),
            axis=1,
        )
        df["p3"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) - 4 + add),
            axis=1,
        )
        df["p2"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) - 3 + add),
            axis=1,
        )
        df["p1"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) - 2 + add),
            axis=1,
        )

        df["p1'"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) - 1 + add),
            axis=1,
        )
        df["p2'"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) + 0 + add),
            axis=1,
        )
        df["p3'"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) + 1 + add),
            axis=1,
        )
        df["p4'"] = df.apply(
            lambda row: self._get_aa(row["Protein"], int(row[start_or_end]) + 2 + add),
            axis=1,
        )
        return df


def calc_height(p, q):
    """
    KL-divergence
    """
    I = (p * np.log2(p / q)).sum(axis=1)
    return (p.T * I).T


def get_height_from_norm(test, back):
    height = calc_height(test.T, back.T)
    return height
