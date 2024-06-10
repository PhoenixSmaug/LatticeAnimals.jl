# LatticeAnimals.jl

This packages generates random lattice animals (or sometimes called [polyforms](https://en.wikipedia.org/wiki/Polyform)) from the standard percolation model using the [Metropolis-Hasting](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm) algorithm.

The polyform is specified by its lattice basis matrix and a list of neighbours, where each entry is a linear combination of basis vectors leading to a neighbouring lattice point. So for example the square tesselation of the polyomino squares has the basis matrix $B$ and neighbours $N$:
```math
B = \begin{pmatrix}1 & 0 \\0 & 1 \end{pmatrix} \quad \quad \quad N = \left\{\begin{pmatrix}1 \\0 \end{pmatrix}, \begin{pmatrix}0 \\1 \end{pmatrix}, \begin{pmatrix}-1 \\0 \end{pmatrix}, \begin{pmatrix}0 \\-1 \end{pmatrix}\right\}
```

The following functions are exported:

* `Poly(n, p, basis, neighbours)`: Generate a random polyform - specified in the basis-neighbours notation - of size n and percolation factor p 

* `Polyomino(n, p)`: Generate a random polyomino of size n and percolation factor p

* `Polyhex(n, p)`: Generate a random polyhex of size n and percolation factor p

* `polyPlot(tiles, basis, path)`: Draw the polyform as a scatter plot and save to path

* `setDimension(d)`: Change the number of dimensions (default value: 2), a Julia restart is required afterwards