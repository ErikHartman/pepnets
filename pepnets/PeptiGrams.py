import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib

plt.rcParams["font.family"] = "Arial"

class PeptiGram:
    def __init__(self, dm, design):
        self.dm = dm
        self.design = design

    def plot_peptigram(
        self,
        protein: str,
        groups: list,
        days: list,
        cmaps: list = None,
        color: dict = None,
        cbar: bool = False,
        print_mods: bool = True,
        xlim: list = [0, 0],
        size_factor: int = 1,
        save_str : str = None
    ):
        data = self.dm
        data["Protein"] = data["Cluster"].apply(lambda x: x.split("_")[0])
        data = data[data["Protein"] == protein].copy().reset_index()
        design = self.design
        design = design[design["group"].isin(groups)]
        design = design[design["day"].isin(days)]

        all_samples = design[design["day"].isin(days)]["sample"].values

        min_start = min(data["Start"].astype(int))
        max_end = max(data["End"].astype(int))

        if xlim == [0, 0]:
            xlim = [min_start - 1, max_end + 1]

        sequence_range = range(min_start, max_end)
        n_bases = len(sequence_range)

        if not cmaps:
            cmaps = [
                "Reds",
                "Purples",
                "Blues",
                "Greens",
                "Oranges",
                "RdPu",
                "YlOrBr",
                "PuRd",
                "OrRd",
                "BuPu",
                "GnBu",
                "PuBu",
                "YlGnBu",
                "PuBuGn",
                "YlOrRd",
                "BuGn",
                "YlGn",
            ]
        cmaps = cmaps * 10

        cluster_idx = sorted(
            data["Cluster"].unique().tolist(), key=lambda x: int(x.split("_")[-1])
        )

        """
        This gets height_ratios
        """
        heights = []
        for group in groups:
            max_heights = []
            samples = design[design["group"] == group]
            for day in days:
                max_height = 20
                for cluster in data["Cluster"].unique():
                    columns = samples[samples["day"] == day][
                        "sample"
                    ].values.tolist() + [
                        "Start",
                        "End",
                        "Cluster",
                    ]
                    datamatrix = data[columns].set_index(["Start", "End", "Cluster"])
                    datamatrix["mean"] = datamatrix.mean(axis=1)
                    datamatrix = datamatrix[datamatrix["mean"] > 0]
                    datamatrix.reset_index(inplace=True)
                    sequence = np.zeros(len(sequence_range))
                    for row in datamatrix.iterrows():
                        row = row[1]
                        start = row["Start"] - min_start
                        end = row["End"] - min_start
                        cluster_id = cluster_idx.index(row["Cluster"])
                        max_i = 0
                        for i in range(int(start), int(end)):
                            sequence[i] += 1
                            max_i = max(max_i, sequence[i])
                            max_height = max(max_i, max_height)
                    max_heights.append(max_height)

            heights.append(max(max_heights))

        if sum(heights) == 0:
            print("No peptides in data")
            return None

        heights.insert(0, sum(heights) / 10)

        fig, axs = plt.subplots(
            len(heights),
            len(days),
            figsize=(
                n_bases * len(days) * size_factor / 45,
                sum(heights[1:]) * size_factor / 18,
            ),
            gridspec_kw={"height_ratios": np.array(heights)},
            sharex=True,
        )

        max_color_int = 0

        max_mean = max(data[all_samples].mean(axis=1))

        for day in days:
            day_index = days.index(day)
            for group in groups:
                group_index = groups.index(group)

                if len(days) == 1:
                    ax = axs[group_index + 1]
                else:
                    ax = axs[group_index + 1, day_index]

                samples = design[design["group"] == group]
                columns = samples[samples["day"] == day]["sample"].values.tolist() + [
                    "Start",
                    "End",
                    "Cluster",
                    "Peptide",
                    "Mods",
                ]

                if len(columns) == 3:
                    sns.despine(top=True, bottom=True, right=True, left=True, ax=ax)
                    ax.set_visible("False")
                    ax.set_yticks([])
                else:
                    datamatrix = (
                        data[columns]
                        .copy()
                        .set_index(["Start", "End", "Cluster", "Peptide", "Mods"])
                    )
                    datamatrix["mean"] = datamatrix.mean(axis=1, numeric_only=True)

                    datamatrix = datamatrix[datamatrix["mean"] > 0].copy()
                    datamatrix.reset_index(inplace=True)

                    spaces = np.zeros((int(max(heights)), int(max_end)))

                    datamatrix["length"] = datamatrix["End"] - datamatrix["Start"]
                    datamatrix.sort_values(
                        ["Start", "length"], inplace=True, ascending=[True, False]
                    )

                    max_height = 0

                    for row in datamatrix.iterrows():
                        row = row[1]
                        start = int(row["Start"])
                        end = int(row["End"])
                        cluster_id = cluster_idx.index(row["Cluster"])
                        cmap = plt.get_cmap(cmaps[cluster_id])

                        if color and isinstance(color, str) and color == "constant":
                            color_intensity = 150
                        elif color:
                            sequence = row["Peptide"]
                            color_intensity = color[sequence]
                        else:
                            color_intensity = 0.25 + (row["mean"] / (max_mean)) / 1.25

                        max_color_int = max(color_intensity, max_color_int)

                        for height in range(spaces.shape[0]):
                            position = spaces[height, start:end]
                            if sum(position) == 0:
                                spaces[height, start : end + 1] = 1
                                ax.plot(
                                    [start + 0.5, end - 0.5],
                                    [-height - 1, -height - 1],
                                    color=cmap(color_intensity),
                                    linewidth=2 * size_factor,
                                    solid_capstyle="round",
                                )
                                if print_mods:
                                    mods = self.parse_mod(row["Mods"].split("|"))
                                    for mod in mods:
                                        mod_index = mod[0]
                                        mod_letter = mod[1]
                                        ax.text(
                                            x=start + mod_index,
                                            y=-height - 1,
                                            s=mod_letter,
                                            fontsize=4,
                                        )
                                break
                    ax.set_ylim([-heights[group_index + 1] - 2, 0])
                    ax.set_xlim(xlim)
                    ax.set_ylabel(group)
                    ax.set_yticks([])
                    ax.axhline(0, 0, 1, color="lightgray", linestyle="--")
                    sns.despine(top=True, bottom=True, right=True, left=True, ax=ax)

            for day in days:
                day_index = days.index(day)
                columns = design[design["day"] == day]["sample"].values.tolist() + [
                    "Start",
                    "End",
                    "Cluster",
                ]
                datamatrix = data[columns].set_index(["Start", "End", "Cluster"]).copy()
                datamatrix["mean"] = datamatrix.mean(axis=1)
                datamatrix.reset_index(inplace=True)

                for cluster in datamatrix["Cluster"].unique():
                    if len(days) == 1:
                        ax = axs[0]
                    else:
                        ax = axs[0, day_index]

                    cluster_id = cluster_idx.index(cluster)
                    cmap = plt.get_cmap(cmaps[cluster_id])
                    spaces = np.zeros(int(max_end))
                    cl = datamatrix[datamatrix["Cluster"] == cluster].copy()
                    for row in cl.iterrows():
                        if row[1]["mean"] > 0:
                            start = int(row[1]["Start"])
                            end = int(row[1]["End"])
                            spaces[start:end] += 1
                    ax.bar(
                        list(range(len(spaces))),
                        height=spaces,
                        color=cmap(0.5),
                        width=1,
                        alpha=0.5,
                    )

                    # ax.annotate(
                    #     cluster.split("_")[-1],
                    #     xy=(np.argmax(spaces), np.max(spaces)),
                    #     size=10,
                    # )
                ax.set_xlim(xlim)
                # ax.set_title(day, pad=5)
                ax.set_yticks([])
                sns.despine(ax=ax, top=True, bottom=True, right=True, left=True)

        plt.suptitle(f"{protein}", y=1 + 0.1 * size_factor, fontsize=20 * size_factor)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.05)

        if cbar:
            cax = fig.add_axes([1.1, 0.1, 0.025, 0.8])
            norm = matplotlib.colors.Normalize(vmin=0, vmax=max_color_int)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            fig.colorbar(sm, cax=cax, label="AMP Score")
        if save_str:
            plt.savefig(save_str, dpi=300)

    def parse_mod(self, mods):
        parsed_mods = []
        mods = [m for m in mods if m != ""]
        for mod in mods:
            mod_list = mod[1:-1].split(",")
            mod_index = int(mod_list[0])
            m = str(mod_list[1])
            parsed_mods.append((mod_index, m))
        return parsed_mods
