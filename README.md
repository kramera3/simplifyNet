# simplifyNet

Package for network sparsification.

## Description

An R package for network sparsification with a variety of novel and known network sparsification techniques. All network sparsification reduce the number of edges, not the number of nodes. A network is usually a large, complex weighted graph obtained from real-world data. It is commonly stored as an adjacency matrix or edge list. Network sparsification is sometimes referred to as network dimensionality reduction.

## Getting Started

Install and load *devtools* package:

``` sh
install.packages("devtools")
```

Use install_github function to pull and install simplifyNet in your session:

``` sh
install_github("kramera3/simplifyNet")
```

### Prerequisites

The following packages are required:

``` sh
igraph, sanic, Matrix, tidyr, methods, fields, stats, dplyr
```

Also set up the working directory:

``` sh
setwd("<em>working directory</em>")
```

## simplifyNet

*simplifyNet* is a **R** package for network sparsification. It contains a suite of different network sparsification algorithms to output a sparsified network.

**Global Network Sparsification:**

Global network sparsification. Uses a threshold cutoff to remove all edges below a certain edge weight or removes a certain proportion of lowest edge weight edges.

``` sh
gns(E_List, remove.prop, cutoff)
```

**Arguments**

-   **E_List:** An edge list of the format \|n1\|n2\|weight\|)
-   **remove.prop:** The proportion of highest weighted edges to retain. A value between 0 and 1.
-   **cutoff:** Threshold value for edge weight thresholding.

**LANS:**

Local Adaptive Network Sparsification from the paper by [Foti et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0016431#s5)

``` sh
lans(Adj, remove.prop, output)
```

**Arguments**

-   **Adj:** Weighted adjacency matrix of the network.
-   **remove.prop:** Alpha value threshold to designate statistically "unimportant" edges by edge weight.
-   **output:** Designates if the output should be directed or undirected. Default is that the output is the same as the input based on adjacency matrix symmetry. If the default value is to be overridden, set as either "undirected" or "directed".

**Sparsification by Edge Effective Resistances:**

Sparsification by sampling edges proportional to their effective resistances as formulated by [Spielman and Srivastava](https://epubs.siam.org/doi/abs/10.1137/080734029?casa_token=2zbhxtOO76wAAAAA:VUhdSEpKiYM2vX3_yEbrOnhSOJaGnXiTjSmlmvmqHP0jb1-sS5tIaF1V5B4UFReBAcRON8WU7Q). This requires two discrete steps: (1) approximating the effective resistances for all edges, (2) sampling them according to the method devised by Spielman and Srivastava.

``` sh
effR = EffR(E_List, epsilon, type="kts", tol)
EffRSparse(n, E_List, q, effR)
```

1.  **EffR**, effective resistances calculator.

    -   **E_List**: Edge list formatted \| n1 \| n2 \| weight \|.

    -   **epsilon**: Governs the relative fidelity of the approximation methods 'spl' and 'kts'. The smaller the value, the greater the fidelity of the approximation and the greater the space and time requirements. Default value is 0.1.

    -   **type**: There are three methods.

        $1$ 'ext' which exactly calculates the effective resistances (WARNING! Not ideal for large graphs).

        $2$ 'spl' which approximates the effective resistances of the inputted graph using the original Spielman-Srivastava algorithm.

        $3$ 'kts' which approximates the effective resistances of the inputted graph using the implementation by Koutis et al. (ideal for large graphs where memory usage is a concern).

    -   **tol**: Tolerance for the linear algebra (conjugate gradient) solver to find the effective resistances. Default value is 1e-10.

2.  **EffRSparse**, network sparsification through sampling effective resistances.

    -   **n**: The number of nodes in the network.

    -   **E_List**: Edge list formatted \| n1 \| n2 \| weight \|.

    -   **q**: The numbers of samples taken. The fidelity to the original network increases as the number of samples increases, but decreases the sparseness.

    -   **effR**: Effective resistances corresponding to each edge. Should be in the same order as "weight".

## Authors

-   **Dr. Andrew Kramer** - Primary Author [simplifyNet](https://github.com/kramera3/simplifyNet)
-   **Alexander Mercier** - Package Maintainer
-   **Shubhankar Tripathi** - Contributor
-   **Tomlin Pulliam** - Contributor
-   **John Drake** - Contributor

## Method Acknowledgements

-   **EffR** and **EffRSparse** are based on work by [Daniel A. Spielman and Nikihl Srivastava](https://epubs.siam.org/doi/abs/10.1137/080734029?casa_token=2zbhxtOO76wAAAAA:VUhdSEpKiYM2vX3_yEbrOnhSOJaGnXiTjSmlmvmqHP0jb1-sS5tIaF1V5B4UFReBAcRON8WU7Q).
-   **EffR** also based on work by [Koutis et al.](https://www.cs.cmu.edu/~jkoutis/papers/stacs239koutis.pdf)
-   **toivonen** based on work by [Hannu Toivonen et al.](https://link.springer.com/chapter/10.1007/978-3-642-13062-5_21)
-   **lans** based on work by [Foti et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3035633/)

## License

-   **GNU General Public License**
