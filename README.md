# Structure Factor

Structure factor is calculated for a given xyz file.

The script reads in input.xyz and inc.params.

## input.xyz
N (row) by 3 (col) in which each row corresponds to a single atom/bead and each col x, y, and z.

## inc.params
* `Lx Ly Lz`, simulation box lengths in x, y, and z axes;
* `n`, number of points in each axis in the reciprocal space;

## example
<img src="https://github.com/suwonbae/structurefactor/blob/main/sq.png" style="width: 350px;"/>
