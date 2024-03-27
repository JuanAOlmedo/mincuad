#include <stdlib.h>
#include <stdio.h>
#include "include/gc.h"

/* Equivalente a init_matrix, solo que la matriz resultante no
 * va a ser guardada para ser liberada más tarde */
static Matrix init_matrix_man(unsigned rows, unsigned cols)
{
    return (Matrix) {rows, cols, malloc(sizeof(double) * rows * cols)};
}

Matrix init_matrix(unsigned rows, unsigned cols)
{
    Matrix result = init_matrix_man(rows, cols);
    push_to_matrix_list(&result);
    
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
    remove_from_matrix_list(M->mat);
    free(M->mat);
    M->mat = NULL;
}

void free_all(void)
{
    double *aux;
    struct matrix_list *matrices = get_matrix_list();

    while (matrices != NULL) {
        aux = matrices->mat;
        remove_from_matrix_list(aux);
        free(aux);
    }
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

    result.f = 0;
    if (A->rows == 1)
        result.f = *A->mat;
    else {
        cof = init_matrix_man(A->rows - 1, A->cols - 1);
        if ((result.error = cof.mat == NULL))
            return result;

        /* Tomar la primera fila y calcular el determinante de
         * los cofactores de cada entrada de la fila */
        for (int i = 0; i < A->rows; i++) {
            cofactor_matrix_of(A, 0, i, &cof);
            det_resu = determinant(&cof);

            if ((result.error = det_resu.error)) {
                free(cof.mat);
                return result;
            }

            if (i % 2 == 0)
                result.f += *(A->mat + i) * det_resu.f;
            else
                result.f -= *(A->mat + i) * det_resu.f;
        }

        free(cof.mat);
    }

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

Matrix transpose_matrix_of(Matrix *A)
{
    Matrix A_t = init_matrix(A->cols, A->rows);
    unsigned i, j;

    if (A_t.mat != NULL)
        for (i = 0; i < A_t.rows; i++)
            for (j = 0; j < A_t.cols; j++)
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
