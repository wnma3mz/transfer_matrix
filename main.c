#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "svd.c"


// Parameters for SVD
#define M 1          // Number of rows in the matrix
#define N 801        // Number of columns in the matrix
#if (M>N)      // MN is the larger number of M and N
#define MN M
#else
#define MN N
#endif

const int num = 801;
const double pi = 3.1416;
const int NN = 10;


double **dmatrix(int, int, int, int);
double *dvector(int, int);
void svdcmp(double **, int, int, double *, double **);


double _Complex **init(int m, int n) {
    /* Initialization of complex matrix
     * m: Number of rows in the matrix
     * n: Number of columns in the matrix
     */
    int i;
    double _Complex **Matrix;
    Matrix = (double _Complex**) malloc(sizeof(double _Complex*) * (m+1));
    for(i = 1; i <= m; i++) {
        Matrix[i] = (double _Complex*) malloc(sizeof(double _Complex) * (n+1));
    }
    return Matrix;
}


void print_r(double **a, int m, int n) {
    /* Print matrix
     * a: Matrix that needs to be output
     * m: Number of rows in the matrix
     * n: Number of columns in the matrix
     */
    int i, j;
    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            printf("%le ", a[i][j]);
        }
        printf("\n");
    }
}

double _Complex **complex_matmul(double _Complex **k1, double _Complex **k2, int r1, int c1, int c2) {
    /* Multiplication of complex matrices
     * k1, k2: Two matrices that need to be multiplied
     * r1: The number of rows of k1
     * c1: The number of columns of k1
     * c2: The number of rows of k2
     */
    int i, j, k;
    complex double sum;
    double _Complex **res;

    res = init(r1, c2);

    for (i = 1; i <= r1; i++) {
        for (j = 1; j <= c2; j++) {
            sum = 0.0 + 0 * _Complex_I;
            for (k = 1; k <= c1; k++) {
                sum += k1[i][k] * k2[k][j];
            }
            res[i][j] = sum;
        }
    }

    return res;
}

double **matmul(double **k1, double **k2, int r1, int c1, int c2) {
    /* Multiplication of complex matrices
     * k1, k2: Two matrices that need to be multiplied
     * r1: The number of rows of k1
     * c1: The number of columns of k1
     * c2: The number of rows of k2
     */
    int i, j, k;
    double sum;
    double **res;

    res = dmatrix(1, MN, 1, MN);
    for (i = 1; i <= r1; i++) {
        for (j = 1; j <= c2; j++) {
            sum = 0.0;
            for (k = 1; k <= c1; k++) {
                sum += k1[i][k] * k2[k][j];
            }
            res[i][j] = sum;
        }
    }
    return res;
}

double **transpose(double **k, int r, int c) {
    /* Transpose of matrix
     * k: Matrix that needs to be transposed
     * r: Number of rows in the matrix
     * c: Number of columns in the matrix
     */
    int i, j;
    double **res;

    res = dmatrix(1, r, 1, c);
    for (i = 1; i <= r; i++) {
        for (j = 1; j <= c; j++) {
            res[i][j] = k[j][i];
        }
    }
    return res;
}

double inv(double **a, double **u, int m, int n, double w[], double **v, double **a2) {
    /* Reference: https://blog.csdn.net/zhuiqiuk/article/details/69390357
     * a: Perform a matrix of SVD
     * u, w, v: Three matrices decomposed by matrix a
     * m, n: Number of rows in the matrix, Number of columns in the matrix
     * a2: Another matrix
     */
    int i,j,k;
    double t;
    double t1[MN], t2[MN];

    for (i = 1; i <= M; i++) {
        for (j = 1; j <= N; j++)
            u[i][j] = a[i][j];
    }

    for (i = 1; i <= M; i++) {
        for (j = 1; j <= N; j++)
            u[i][j] = a[i][j];
    }

    svdcmp(u, MN, MN, w, v);
    /* Sort the singular values in descending order */
    for (i = 1; i <= N; i++) {
        for (j = i+1; j <= N; j++) {
            if (w[i] < w[j]) {
                t = w[i];
                w[i] = w[j];
                w[j] = t;
                /* Column position of the exchange matrix U, V */
                for (k = 1; k <= M; k++) {
                    t1[k] = u[k][i];
                }
                for (k = 1; k <= M; k++) {
                    u[k][i] = u[k][j];
                }
                for (k = 1; k <= M; k++) {
                    u[k][j] = t1[k];
                }

                for (k = 1; k <= N; k++) {
                    t2[k] = v[k][i];
                }
                for (k = 1; k <= N; k++) {
                    v[k][i] = v[k][j];
                }
                for (k = 1; k <= N; k++) {
                    v[k][j] = t2[k];
                }
            }
        }
    }


    /* U: MxM */
    /* M singular values */
    /* V: NxN */

    double **W;
    W = dmatrix(1, MN, 1, MN);

    for (i = 1; i<= M; i++) {
        for (j =1; j <= N; j++) {
            if (i==j) {
                W[i][j] = w[i];
            } else {
                W[i][j] = 0.0;
            }
        }
    }


    /* 计算B=V*(1/W)*U.T */
    // B is the result matrix
    double **temp, **V, **B, **U;
    V = dmatrix(1, MN, 1, MN);
    U = dmatrix(1, MN, 1, MN);
    temp = dmatrix(1, MN, 1, MN);
    B = dmatrix(1, MN, 1, MN);

    double sum;

    V = v;
    U = transpose(u, M, M);

    // Since there is only one singular value, turn it into its own reciprocal
    W[1][1] = 1 / W[1][1];


    temp = matmul(V, W, N, N, M);
    B = matmul(temp, U, N, M, M);

    // Calculate a2/a. Matrix multiplication, ie vector multiplication
    double res = 0.0;
    for (i = 1; i <= N; i++) {
        res += B[i][1] * a2[1][i];
    }

    return res;
}

double k_div(double **k1, double **k2) {

    double *w;
    double **u, **v;

    u = dmatrix(1, MN, 1, MN);
    w = dvector(1, MN);
    v = dmatrix(1, MN, 1, MN);

    return inv(k1, u, MN, MN, w, v, k2);
}
int main() {

    int i, j, k;


    double n1=2.33;
    double n2=1.50;	    // 两层介质的折射率
    double h1=200e-9;
    double h2=100e-9;	// 介质2厚度为0
    int A=1;
    int B=0;


    double _Complex **M1, **M2;
    double T[num];
    double **k1, **k2;

    int wl[num];
    double e1[num], e2[num];

    complex double c1, c2;

    for (i = 0; i < num; i++) {
        wl[i] = i + 100;
    }

    k1 = dmatrix(1, MN, 1, MN);
    k2 = dmatrix(1, MN, 1, MN);

    for (i = 0; i < num; i++) {
        k1[1][i+1] = 2*pi*n1 / (wl[i]*1e-9);
        e1[i] = 2*pi*n1*h1 / (wl[i]*1e-9);
        k2[1][i+1] = 2*pi*n2 / (wl[i]*1e-9);
        e2[i] = 2*pi*n2*h2 / (wl[i]*1e-9);
    }

    // Calculate k1/k2 and k2/k1
    double temp1, temp2;
    temp1 = k_div(k2, k1);
    temp2 = k_div(k1, k2);
    c1 = 0.5 * _Complex_I * (temp1 + temp2);
    c2 = 0.5 * _Complex_I * (temp1 - temp2);

    printf("k1/k2:%e\n", temp1);
    printf("k2/k1:%e\n", temp2);
    printf("c1:%.3f+%.3ei\n", creal(c1), cimag(c1));
    printf("c2:%.3f+%.3ei\n", creal(c2), cimag(c2));

    M1 = init(2,2);
    M2 = init(2,2);

    for (i = 0; i < num; i++) {
        M1[1][1] = cexp( _Complex_I * e1[i]) * (cos(e2[i]) + c1 * sin(e2[i]));
        M1[1][2] = cexp(-_Complex_I * e1[i]) * (-c2 * sin(e2[i]));
        M1[2][1] = cexp( _Complex_I * e1[i]) * c2 * sin(e2[i]);
        M1[2][2] = cexp(-_Complex_I * e1[i]) * (cos(e2[i]) - c1 * sin(e2[i]));
        M2 = M1;
        // Matrix multiply 10-1 times
        for (k = 1; k < NN; k++) {
            M2 = complex_matmul(M2, M1, 2, 2, 2);
        }
        T[i] = pow(cabs((M2[1][1] * M2[2][2] - M2[2][1] * M2[1][2]) / ((A * M2[2][2] - B * M2[1][2]))), 2);
    }

    // Output the first 100 values of T for inspection
    for (i = 0; i < 100; i++) {
        printf("%.3e\t", T[i]);
    }

    return 0;
}


