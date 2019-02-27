## Transfer_matrix

Received help from kabume, re-implementing the **Transfer_Matrix** matlab code in C language.

Original project address: https://github.com/kabume/Transfer_matrix_1D_PHC

Highlight content:

1. Achieved pseudo-inverse(Generalized inverse)

2. Using *complex.h* and multiply the complex matrix

### main.m

Implemented **Transfer_Matrix** using matlab code.

### svd.c

The *svd.h* is used to complete the SVD matrix decomposition, and then the pseudo-inverse matrix is solved. From http://cacs.usc.edu/education/phys516/src/TB/svdcmp.c.

### main.c

Implement the main program of **Transfer_Matrix**.

### Usage

**Environment**: gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-28)

```bash
gcc main.c -lm -o main
```

