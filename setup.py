from setuptools import setup

requirements=[
    "networkx",
    "leidenalg",
    "pandas",
    "matplotlib",
    "numpy",
    "seaborn"
    "igraph",
    "logomaker"
]

setup(
    author="Erik Hartman",
    author_email="erik.hartman@med.lu.se",
    name="pepnets",
    version="0.0.1",
    packages=["pepnets"],
    requires=requirements
)
