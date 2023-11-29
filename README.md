<p align="center">
    <img src="logo.png", width="250" />
<p>

# pepnets: clustering peptides on protein backbones

This package clusters linear peptides that are degradation products of proteins based on their sequence.

## install
To install the package, clone this repo:
```
> git clone https://github.com/ErikHartman/pepnets.git
```
and then pip install:
```
> pip install path_to_root
```
## usage


```py
pnet = PeptideNetwork(
    datamatrix=datamatrix,
    protein_database=sus_scrofa,
)

pnet.get_clusters(resolution=0.8)
```
```py
pnet.plot_protein(protein="APOA1")
```
![network](plots/APOA1.png "network")
```py

fm = FeatureMatrix(datamatrix, pnet.clusters)
```
```py
design = pd.read_csv("../data/design.csv")
peptigram = PeptiGram(fm.datamatrix, design)
peptigram.plot_peptigram(
    "APOA1", groups=["S. aureus", "P. aeruginosa"], days=["Day 1"], size_factor=0.3,
)
```
![peptigram](plots/APOA1_pg.png "peptigram")

For more examples, see `notebooks`.

## dependencies
This package is heavily dependent on `networkx` for graph creation, `leidenalg` for community detection, and several packages for manipulation and plotting (`pandas, numpy, matplotlib, seaborn`).