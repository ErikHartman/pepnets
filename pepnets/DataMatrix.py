import pandas as pd
import numpy as np


class DataMatrix:
    def __init__(self, data: pd.DataFrame, design: pd.DataFrame):
        self.data = data.replace(0, np.nan)
        self.design = design
        self.groups = self.design["group"].unique().tolist()
        self.X = data.values

    def get_group(self, group: str):
        if group not in self.design["group"].values.tolist():
            print(f"{group} is not in data. {self.design['group'].unique().tolist()}")
        samples = self.design[self.design["group"] == group]["sample"].values.tolist()
        sample_data = self.data[samples]

        return sample_data

    def get_missing_value_fractions(self):
        all_n_vals = len(self.data.index) * len(self.data.columns)
        all_na_vals = self.data.isna().sum().sum()
        fractions = {}
        for group in self.groups:
            group_data = self.get_group(group)
            group_data = group_data.dropna(axis=0, how="all")
            group_data.replace(0, np.nan)
            total_n_vals = len(group_data.index) * len(group_data.columns)
            na_vals = group_data.isna().sum().sum()
            fractions[group] = na_vals / total_n_vals
        fractions["All"] = all_na_vals / all_n_vals
        return fractions

    def get_mean_intensities(self):
        df_dict = {}
        for group in self.groups:
            group_data = self.get_group(group).replace(0, np.nan)
            means = group_data.mean(axis=1)
            df_dict[group] = means
        return pd.DataFrame(df_dict)

    def get_flat_data(self):
        flat_matrix = pd.DataFrame()
        for group in self.groups:
            group_data = self.get_group(group)
            group_data = group_data.replace(0, np.nan).stack().reset_index()
            print(group_data)
