# Structure Factor

Written in Fortran.

## input.xyz
N (row) by 3 (col) in which each row corresponds to a single atom/bead and each col x, y, and z.

## inc.parameters
* `num_atoms`, number of total atoms/beads = N;
* `n`, number of points in each axis in the reciprocal space;
* `Lx`, simulation box length in x axis;
* `Ly`, length in y axis;
* `Lz`, length in z axis;
