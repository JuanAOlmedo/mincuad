#include "include/matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define pointer_to(A, row, col) ((A).mat + (row) * (A).cols + col)

Matrix matrix_new(unsigned rows, unsigned cols)
{
    Matrix result = (Matrix) {rows, cols, malloc(sizeof(double) * rows * cols)};

    return result;
}

int matrix_is_valid(Matrix A)
{
    return A.mat == NULL;
}

Matrix matrix_copy(Matrix A)
{
    Matrix B = matrix_new(A.rows, A.cols);

    memcpy(B.mat, A.mat, A.cols * A.rows * sizeof(double));

    return B;
}

double matrix_read(Matrix A, unsigned row, unsigned col)
{
    if (row < A.rows && col < A.cols)
        return *pointer_to(A, row, col);

    return 0;
}

void matrix_write(Matrix *A, unsigned row, unsigned col, float x)
{
    if (row < A->rows && col < A->cols)
        *pointer_to(*A, row, col) = x;
}

double *matrix_col(Matrix A, unsigned col)
{
    if (col < A.cols) {
        double *col_arr = malloc(sizeof(double) * A.rows);
        for (unsigned i = 0; i < A.rows; i++)
            col_arr[i] = *pointer_to(A, i, col);

        return col_arr;
    }

    return NULL;
}

double *matrix_row(Matrix A, unsigned row)
{
    if (row < A.rows) {
        double *row_arr = malloc(sizeof(double) * A.cols);
        memcpy(row_arr, pointer_to(A, row, 0), sizeof(double) * A.cols);

        return row_arr;
    }

    return NULL;
}

void matrix_free(Matrix *A)
{
    free(A->mat);
    A->mat = NULL;
}

Matrix matrix_smultiply(Matrix A, double x)
{
    Matrix result = matrix_new(A.rows, A.cols);

    for (unsigned i = 0; i < result.rows * result.cols; i++)
        result.mat[i]= A.mat[i] * x;

    return result;
}

Matrix matrix_multiply(Matrix A, Matrix B)
{
    Matrix result = matrix_new(A.rows, B.cols);
    unsigned i, j, k;
    double *result_i_j = result.mat;

    if (A.cols != B.rows) {
        if (result.mat != NULL)
            matrix_free(&result);
        return result;
    }

    for (i = 0; i < result.rows; i++)
        for (j = 0; j < result.cols; j++, result_i_j++)
            for (*result_i_j = k = 0; k < A.cols; k++)
                *result_i_j += *pointer_to(A, i, k) * *pointer_to(B, k, j);

    return result;
}

Matrix matrix_add(Matrix A, Matrix B)
{
    Matrix result = (Matrix) {0, 0, NULL};

    if (A.rows == B.rows && A.cols == B.cols) {
        result = matrix_new(A.rows, A.cols);
        for (unsigned i = 0; i < result.rows * result.cols; i++)
            result.mat[i]= A.mat[i] + B.mat[i];
    }

    return result;
}

void matrix_print(Matrix A)
{
    unsigned i, j;

    for (i = 0; i < A.rows; i++) {
        for (j = 0; j < A.cols; j++)
            printf("%10.5f ", *pointer_to(A, i, j));

        putchar('\n');
    }
}

void matrix_get(Matrix *A)
{
    for (double *mat = A->mat; (mat - A->mat) < (A->rows * A->cols); mat++)
        scanf("%lf", mat);
}

/* Devuelve el signo de la permutación p */
static int sgnp(unsigned *p, unsigned n)
{
    int s = 1;
    for (unsigned i = 0; i < n; i++) {
            if (p[i] % 2 == 0)
                s *= - 1;

            for (unsigned j = i + 1; j < n; j++)
                if (p[j] > p[i])
                    p[j] -= 1;
    }
    return s;
}

double matrix_det(Matrix A)
{
    if (A.rows != A.cols)
        return 0;

    struct LU decomp = matrix_lu(A);
    double result = sgnp(decomp.p, decomp.n);

    for (unsigned i = 0; i < decomp.n; i++)
        result *= *pointer_to(decomp.U, i, i);

    matrix_free(&decomp.L);
    matrix_free(&decomp.U);
    free(decomp.p);
    return result;
}

static void forward_sub(Matrix L, Matrix b, Matrix *x, unsigned *p)
{
    x->mat[0] = *pointer_to(b, p[0], 0) / *pointer_to(L, 0, 0);
    for (unsigned k = 1; k < x->rows; k++) {
        x->mat[k] = 0;
        for (unsigned j = 0; j < k; j++)
            x->mat[k] += *pointer_to(*x, j, 0) * *pointer_to(L, k, j);

        x->mat[k] = (b.mat[p[k]] - x->mat[k]) / *pointer_to(L, k, k);
    }
}

static void back_sub(Matrix U, Matrix b, Matrix *x)
{
    unsigned n = U.rows;
    x->mat[n - 1] = *pointer_to(b, n - 1, 0) / *pointer_to(U, n - 1, n - 1);
    for (int k = n - 2; k >= 0; k--) {
        x->mat[k] = 0;
        for (unsigned j = k + 1; j < n; j++)
            x->mat[k] += *pointer_to(*x, j, 0) * *pointer_to(U, k, j);

        x->mat[k] = (b.mat[k] - x->mat[k]) / *pointer_to(U, k, k);
    }
}

Matrix matrix_inv(Matrix A)
{
    if (matrix_det(A) == 0)
        return (Matrix) {0, 0, NULL};

    Matrix result = matrix_new(A.rows, A.cols),
           b = matrix_new(A.rows, 1),
           w = matrix_new(A.rows, 1),
           x = matrix_new(A.rows, 1);
    unsigned i, j, n = A.rows;

    struct LU decomp = matrix_lu(A);
    for (i = 0; i < n; i++)
        b.mat[i] = 0;

    for (i = 0; i < n; i++) {
        /* b es la columna i-ésima de la identidad */
        b.mat[i] = 1;
        if (i != 0)
            b.mat[i - 1] = 0;

        forward_sub(decomp.L, b, &w, decomp.p);
        back_sub(decomp.U, w, &x);

        for (j = 0; j < n; j++)
            *pointer_to(result, j, i) = x.mat[j];
    }
    matrix_free(&decomp.L);
    matrix_free(&decomp.U);
    free(decomp.p);
    matrix_free(&b);
    matrix_free(&w);
    matrix_free(&x);

    return result;
}

static void swap(unsigned *p, unsigned i, unsigned j)
{
    unsigned aux = p[i];
    p[i] = p[j];
    p[j] = aux;
}

static void permute_rows(Matrix *A, unsigned *p)
{
    unsigned n = A->rows, m = A->cols;
    unsigned q[n];
    memcpy(q, p, sizeof(unsigned) * n);

    double aux[m];
    for (unsigned i = 0; i < n; i++) {
        if (q[i] != i) {
            memcpy(aux, pointer_to(*A, i, 0), sizeof(double) * m);
            memcpy(pointer_to(*A, i, 0), pointer_to(*A, p[i], 0), sizeof(double) * m);
            memcpy(pointer_to(*A, p[i], 0), aux, sizeof(double) * m);
            swap(q, i, q[i]);
        }
    }
}

Matrix matrix_system_solve(Matrix A, Matrix b)
{
    if (A.rows != b.rows || A.rows != A.cols || b.cols != 1)
        return (Matrix) {0, 0, NULL};

    struct LU decomp = matrix_lu(A);

    Matrix x = matrix_new(b.rows, 1),
           w = matrix_new(b.rows, 1);
    forward_sub(decomp.L, b, &w, decomp.p);
    back_sub(decomp.U, w, &x);

    matrix_free(&decomp.L);
    matrix_free(&decomp.U);
    free(decomp.p);
    matrix_free(&w);

    return x;
}

Matrix matrix_transpose_of(Matrix A)
{
    Matrix A_t = matrix_new(A.cols, A.rows);

    for (unsigned i = 0; i < A_t.rows; i++)
        for (unsigned j = 0; j < A_t.cols; j++)
            *pointer_to(A_t, i, j) = *pointer_to(A, j, i);

    return A_t;
}

void matrix_transpose(Matrix *A)
{
    unsigned i, j;
    double aux;

    for (i = 0; i < A->rows; i++)
        for (j = i + 1; j < A->cols; j++) {
            aux = *pointer_to(*A, i, j);
            *pointer_to(*A, i, j) = *pointer_to(*A, j, i);
            *pointer_to(*A, j, i) = aux;
        }
}

static unsigned max_col(Matrix A, unsigned *p, unsigned col)
{
    double max = fabs(*pointer_to(A, p[col], col));
    unsigned max_i = col;

    for (unsigned i = col + 1; i < A.rows; i++)
        if (fabs(*pointer_to(A, p[i], col)) > max) {
            max = fabs(*pointer_to(A, p[i], col));
            max_i = i;
        }

    return max_i;
}

static struct LU lu_man(Matrix A)
{
    Matrix U = matrix_copy(A);
    unsigned n = A.rows,
             *p = malloc(sizeof(unsigned) * n),
             *p_inv = malloc(sizeof(unsigned) * n);

    for (unsigned i = 0; i < n; i++)
        p[i] = p_inv[i] = i;

    for (unsigned k = 0; k < n - 1; k++) {
        unsigned m = max_col(U, p, k);

        double diag = *pointer_to(U, p[m], k);
        if (diag != 0) {
            swap(p_inv, p[m], p[k]);
            swap(p, m, k);
            diag = *pointer_to(U, p[k], k);

            for (unsigned i = k + 1; i < n; i++) {
                double *multiplier = pointer_to(U, p[i], k);
                *multiplier /= diag;
                for (unsigned j = k + 1; j < n; j++)
                    *pointer_to(U, p[i], j) -= *multiplier * *pointer_to(U, p[k], j);
            }
        }
    }
    permute_rows(&U, p_inv);
    Matrix L = matrix_copy(U);
    for (unsigned i = 0; i < n; i++) {
        *pointer_to(L, i, i) = 1;
        for (unsigned j = i + 1; j < n; j++) {
            *pointer_to(L, i, j) = 0;
            *pointer_to(U, j, i) = 0;
        }
    }
    free(p_inv);

    return (struct LU) {L, U, p, n};
}
