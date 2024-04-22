import pandas as pd


def get_database(database: pd.DataFrame, key_col="Entry Name", value_col="Sequence"):
    database_dict = {}
    for _, row in database.iterrows():
        database_dict[row[key_col]] = row[value_col]
    return database_dict


def custom_mode(data, equal_strategy="min"):
    """
    Calculate the mode of a list of integers with a specified strategy for handling equal frequencies.

    Parameters:
        data (list of int): The list of integers.
        equal_strategy (str, optional): Strategy for handling equal frequencies. Possible values are "min" (default) and "max". If "min", the smallest value among those with the maximum frequency is returned. If "max", the largest value among those with the maximum frequency is returned.
    """
    counts = {}
    for value in data:
        counts[value] = counts.get(value, 0) + 1

    max_count = max(counts.values())
    modes = [key for key, count in counts.items() if count == max_count]

    if equal_strategy == "min":

        return min(modes)

    elif equal_strategy == "max":

        return max(modes)
