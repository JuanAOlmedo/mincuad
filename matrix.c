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

Matrix multiply_matrix(Matrix *A, Matrix *B)
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

Matrix copy_matrix(Matrix A)
{
    Matrix B = init_matrix(A.cols, A.rows);

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
                *read_matrix_at(&A_t, i, j) = *read_matrix_at(A, j, i);

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
