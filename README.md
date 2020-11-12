# SimplifyNet

Script for network sparsification

## Getting Started

Install and load _devtools_ package:
```sh
install.packages("devtools")
```

Use install_github function to pull and install simplifyNet in your session:
```sh
install_github("shubhankar/simplifyNet")
```
<strong>OR</strong>
```sh
install_github("kramer3/simplifyNet")
```

### Prerequisites

The following packages are required:
```sh
igraph, cPCG, Matrix, tidyr, fields
```

Also set up the working directory:
```sh
setwd("<em>working directory</em>")
```

## simplifyNet
```
simplifyNet(data, method="Toivonen", model, func, cutoff, remove.prop, num.samples, epsilon, matrix.sparse = FALSE, num.nodes = NULL)
```

**Arguments**

* **data:** Inputted network - can be in format adjacency matrix (matrix or sparseMatrix) or edge list (data.frame |n1|n2|weight|)
* **method:** Either 'Toivonen', 'GlobalSparse', 'LocalAdapt', or 'EffectiveResistance'
* **model:** For iterative refitting, model to evaluate edge importance
* **func:** For iterative refitting, function to score edge importance
* **cutoff:** For 'GlobalSparse', global minimum edge weight threshold
* **remove.prop:** For 'LocalAdapt', proportion of edges to remove from inputted network
* **num.samples:** For 'EffectiveResistance', number of samples to take from edge probability distribution
* **epsilon:** For 'EffectiveResistance', degree to approximate effective resistances
* **matrix.sparse = FALSE:** Return sparseMatrix class object
* **num.nodes = NULL:** If number of nodes must be specified

## Authors

* **Dr. Andrew Kramer** - [simplifyNet](https://github.com/kramer3/simplifyNet)

## License

* **GNU General Public License**


