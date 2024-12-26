#include "include/matrix.h"
#include "include/gc.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Equivalente a init_matrix, solo que la matriz resultante no
 * va a ser guardada para ser liberada más tarde */
static Matrix init_matrix_man(unsigned rows, unsigned cols)
{
    return (Matrix) {rows, cols, malloc(sizeof(double) * rows * cols)};
}

Matrix init_matrix(unsigned rows, unsigned cols)
{
    Matrix result = init_matrix_man(rows, cols);
    GC_push(result.mat, free);

    return result;
}

double *read_matrix_at(Matrix *M, unsigned row, unsigned col)
{
    if (row < M->rows && col < M->cols)
        return M->mat + row * M->cols + col;

    return NULL;
}

double *col(Matrix *M, unsigned col)
{
    if (col < M->cols)
        return M->mat + col;

    return NULL;
}

double *row(Matrix *M, unsigned row)
{
    if (row < M->rows)
        return M->mat + row * M->cols;

    return NULL;
}

void free_matrix(Matrix *M)
{
    GC_remove(M->mat);
    M->mat = NULL;
}

Matrix multiply_matrix(Matrix *A, double x)
{
    Matrix result = init_matrix(A->rows, A->cols);

    if (result.mat != NULL)
        for (unsigned i = 0; i < result.rows * result.cols; i++)
            result.mat[i]= A->mat[i] * x;

    return result;
}

Matrix multiply_matrices(Matrix *A, Matrix *B)
{
    Matrix result = init_matrix(A->rows, B->cols);
    unsigned i, j, k;
    double *result_i_j = result.mat;

    if (result.mat == NULL || A->cols != B->rows) {
        if (result.mat != NULL)
            free_matrix(&result);
        return result;
    }

    for (i = 0; i < result.rows; i++)
        for (j = 0; j < result.cols; j++, result_i_j++)
            for (*result_i_j = k = 0; k < A->cols; k++)
                *result_i_j += *read_matrix_at(A, i, k) * *read_matrix_at(B, k, j);

    return result;
}

Matrix add_matrices(Matrix *A, Matrix *B)
{
    Matrix result = (Matrix) {0, 0, NULL};

    if (A->rows == B->rows && A->cols == B->cols)
        result = init_matrix(A->rows, A->cols);

    if (result.mat != NULL)
        for (unsigned i = 0; i < result.rows * result.cols; i++)
            result.mat[i]= A->mat[i] + B->mat[i];

    return result;
}

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

/* Devuelve en B la Matriz A sin la columna col ni la fila row.
 * Se asume que B está bien inicializada */
static void cofactor_matrix_of(Matrix *A, unsigned row, unsigned col, Matrix *B)
{
    double *matB = B->mat;
    double *matA = A->mat;

    while (matB - B->mat < B->rows * B->cols)
        /* Si matA apunta a la fila row o a la columna col */
        if ((matA - A->mat) % A->cols == col || (matA - A->mat) / A->cols == row)
            matA++; /* No copiar la entrada a B */
        else
            *(matB++) = *(matA++); /* Si no, copiar la entrada a B */
}

ResuF determinant(Matrix *A)
{
    ResuF result, det_resu;
    Matrix cof;

    if ((result.error = A->rows != A->cols))
        return result;

    struct LU decomp = lu(*A);
    result.f = sgnp(decomp.p, decomp.n);

    for (unsigned i = 0; i < decomp.n; i++)
        result.f *= *read_matrix_at(&decomp.U, i, i);

    return result;
}

Matrix invert_matrix(Matrix *A)
{
    ResuF detA = determinant(A), detTemp;
    Matrix result = init_matrix(A->rows, A->cols),
           tempMat = init_matrix_man(A->rows - 1, A->cols - 1);
    unsigned i, j;

    /* Return if det(A) == 0 */
    if (detA.error || detA.f == 0 || result.mat == NULL || tempMat.mat == NULL) {
        if (tempMat.mat != NULL)
            free(tempMat.mat);
        if (result.mat != NULL)
            free_matrix(&result);

        return result;
    }

    for (i = 0; i < A->rows; i++)
        for (j = 0; j < A->cols; j++) {
            cofactor_matrix_of(A, j, i, &tempMat);
            detTemp = determinant(&tempMat);

            if (detTemp.error) {
                free_matrix(&result);
                free(tempMat.mat);
                return result;
            }

            *read_matrix_at(&result, i, j) =
                (((i + j) % 2 == 0) ? 1 : -1) * detTemp.f / detA.f;
        }

    free(tempMat.mat);
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
            memcpy(aux, row(A, i), sizeof(double) * m);
            memcpy(row(A, i), row(A, p[i]), sizeof(double) * m);
            memcpy(row(A, p[i]), aux, sizeof(double) * m);
            swap(q, i, q[i]);
        }
    }
}

static Matrix forward_sub(Matrix L, Matrix b)
{
    Matrix x = init_matrix(L.rows, 1);
    *read_matrix_at(&x, 0, 0) = *read_matrix_at(&b, 0, 0) / *read_matrix_at(&L, 0, 0);
    for (unsigned k = 1; k < x.rows; k++) {
        *read_matrix_at(&x, k, 0) = 0;
        for (unsigned j = 0; j < k; j++)
            *read_matrix_at(&x, k, 0) += *read_matrix_at(&x, j, 0) * *read_matrix_at(&L, k, j);

        *read_matrix_at(&x, k, 0) = (*read_matrix_at(&b, k, 0) - *read_matrix_at(&x, k, 0)) / *read_matrix_at(&L, k, k);
    }
    return x;
}

static Matrix back_sub(Matrix U, Matrix b)
{
    unsigned n = U.rows;
    Matrix x = init_matrix(n, 1);
    *read_matrix_at(&x, n - 1, 0) = *read_matrix_at(&b, n - 1, 0) / *read_matrix_at(&U, n - 1, n - 1);
    for (int k = n - 2; k >= 0; k--) {
        *read_matrix_at(&x, k, 0) = 0;
        for (unsigned j = k + 1; j < n; j++)
            *read_matrix_at(&x, k, 0) += *read_matrix_at(&x, j, 0) * *read_matrix_at(&U, k, j);

        *read_matrix_at(&x, k, 0) = (*read_matrix_at(&b, k, 0) - *read_matrix_at(&x, k, 0)) / *read_matrix_at(&U, k, k);
    }
    return x;
}

Matrix solve(Matrix A, Matrix b)
{
    b = copy_matrix(b);
    unsigned i, j;

    if (A.rows != b.rows || A.rows != A.cols || b.cols != 1)
        return (Matrix) {0, 0, NULL};

    struct LU decomp = lu(A);
    permute_rows(&b, decomp.p);
    Matrix x = forward_sub(decomp.L, b);
    free_matrix(&b);
    b = x;
    x = back_sub(decomp.U, b);
    free_matrix(&b);
    free_matrix(&decomp.L);
    free_matrix(&decomp.U);
    GC_remove(decomp.p);

    return x;
}

Matrix copy_matrix(Matrix A)
{
    Matrix B = init_matrix(A.rows, A.cols);

    if (B.mat != NULL)
        memcpy(B.mat, A.mat, A.cols * A.rows * sizeof(double));

    return B;
}

Matrix transpose_matrix_of(Matrix A)
{
    Matrix A_t = init_matrix(A.cols, A.rows);

    if (A_t.mat != NULL)
        for (unsigned i = 0; i < A_t.rows; i++)
            for (unsigned j = 0; j < A_t.cols; j++)
                *read_matrix_at(&A_t, i, j) = *read_matrix_at(&A, j, i);

    return A_t;
}

void transpose_matrix(Matrix *A)
{
    unsigned i, j;
    double aux;

    for (i = 0; i < A->rows; i++)
        for (j = i + 1; j < A->cols; j++) {
            aux = *read_matrix_at(A, i, j);
            *read_matrix_at(A, i, j) = *read_matrix_at(A, j, i);
            *read_matrix_at(A, j, i) = aux;
        }
}

void print_matrix(Matrix *A)
{
    unsigned i, j;

    for (i = 0; i < A->rows; i++) {
        for (j = 0; j < A->cols; j++)
            printf("%10.5f ", *read_matrix_at(A, i, j));

        putchar('\n');
    }
}

void get_matrix(Matrix *A)
{
    double *mat;

    /* Scan rows * cols doubles from stdin */
    for (mat = A->mat; (mat - A->mat) < (A->rows * A->cols); mat++)
        scanf("%lf", mat);
}

static unsigned max_col(Matrix A, unsigned *p, unsigned col)
{
    double max = fabs(*read_matrix_at(&A, p[col], col));
    unsigned max_i = col;

    for (unsigned i = col + 1; i < A.rows; i++)
        if (fabs(*read_matrix_at(&A, p[i], col)) > max) {
            max = fabs(*read_matrix_at(&A, p[i], col));
            max_i = i;
        }

    return max_i;
}

struct LU lu(Matrix A)
{
    Matrix U = copy_matrix(A);
    unsigned n = A.rows;

    unsigned *p = malloc(sizeof(unsigned) * n);
    GC_push(p, free);

    for (unsigned i = 0; i < n; i++)
        p[i] = i;

    for (unsigned k = 0; k < n - 1; k++) {
        swap(p, max_col(U, p, k), k);

        double diag = *read_matrix_at(&U, p[k], k);
        if (diag != 0) {
            for (unsigned i = k + 1; i < n; i++) {
                double *multiplier = read_matrix_at(&U, p[i], k);
                *multiplier /= diag;
                for (unsigned j = k + 1; j < n; j++) {
                    *read_matrix_at(&U, p[i], j) -= *multiplier * *read_matrix_at(&U, p[k], j);
                }
            }
        }
    }
    permute_rows(&U, p);
    Matrix L = copy_matrix(U);
    for (unsigned i = 0; i < n; i++) {
        *read_matrix_at(&L, i, i) = 1;
        for (unsigned j = i + 1; j < n; j++) {
            *read_matrix_at(&L, i, j) = 0;
            *read_matrix_at(&U, j, i) = 0;
        }
    }
    return (struct LU) {L, U, p, n};
}
