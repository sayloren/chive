## UCSF BMI 203 Algorithms Homework 2018

[![Build
Status](https://travis-ci.org/sayloren/chive.svg?branch=master)](https://travis-ci.org/sayloren/chive)

[Travis Build Results](https://travis-ci.org/sayloren/chive)

:see_no_evil: :hear_no_evil: :speak_no_evil:

#### Clustering

## project file structure

The main file that you will need to modify is `cluster.py` and the corresponding `test_cluster.py`. `utils.py` contains helpful classes that you can use to represent Active Sites. `io.py` contains some reading and writing files for interacting with PDB files and writing out cluster info.

```
.
├── README.md
├── data
│   ...
├── hw2skeleton
│   ├── __init__.py
│   ├── __main__.py
│   ├── cluster.py
│   ├── io.py
│   └── utils.py
└── test
    ├── test_cluster.py
    └── test_io.py
```

## usage

To use the package, first run

```
conda install --yes --file requirements.txt
```

to install all the dependencies in `requirements.txt`. Then the package's
main function (located in `hw2skeleton/__main__.py`) can be run as
follows

```
python -m hw2skeleton -P data test.txt
```

## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.

![a](/images/Cluster_hist.png)

| Cluster | k | silhouette score| rand cluster | rand index |
| ---------- |:----------:|:----------:|:----------:|:----------:|
| Partitioning | 5 | 0.265 | Hierarchical | 0.373 |
| Hierarchical | 5 | 0.407 | Random | 0.008 |
| Random | 5 | -0.125 | partitioning | -0.001 |

## contributors

Original design by Scott Pegg. Refactored and updated by Tamas Nagy.
