#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

static struct matrix_array {
    struct matrix_array *nxt;
    float *mat;
} *matrices = NULL;

static int push_to_matrix_list(Matrix *A);

static void remove_from_matrix_list(float *mat);

ResuM init_matrix(unsigned rows, unsigned cols)
{
    ResuM result;
    
    result.m.rows = rows;
    result.m.cols = cols;
    result.m.mat = malloc(sizeof(float) * rows * cols);
    result.error = result.m.mat == NULL || !push_to_matrix_list(&result.m);
    
    return result;
}

float *read_matrix_at(Matrix *M, unsigned row, unsigned col)
{
    if (row < M->rows && col < M->cols)
        return M->mat + row * M->cols + col;

    return NULL;
}

float *col(Matrix *M, unsigned col)
{
    if (col < M->cols)
        return M->mat + col;

    return NULL;
}

float *row(Matrix *M, unsigned row)
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
    float *aux;

    while (matrices != NULL) {
        aux = matrices->mat;
        remove_from_matrix_list(aux);
        free(aux);
    }
}

static void remove_from_matrix_list(float *mat)
{
    struct matrix_array *aux, *it = matrices;

    if (it->mat != mat) {
        for (; it->nxt != NULL && it->nxt->mat != mat; it = it->nxt)
            ;

        aux = it->nxt;
        it->nxt = it->nxt->nxt;
    } else {
        aux = matrices;
        matrices = matrices->nxt;
    }

    free(aux);
}

static int push_to_matrix_list(Matrix *A)
{
    struct matrix_array *it, *el = malloc(sizeof(struct matrix_array));    

    if (el == NULL)
        return 0;
    
    if (matrices == NULL)
        matrices = el;
    else {
        for (it = matrices; it->nxt != NULL; it = it->nxt)
            ;
        it->nxt = el;
    }

    el->mat = A->mat;
    el->nxt = NULL;

    return 1;
}

ResuM multiply_matrix(Matrix *A, Matrix *B)
{
    ResuM result = init_matrix(A->rows, B->cols);
    int i, j, k;
    float *result_i_j = result.m.mat;

    if ((result.error = result.error || A->cols != B->rows)) {
        if (result.m.mat != NULL)
            free_matrix(&result.m);
        return result;
    }

    for (i = 0; i < result.m.rows; i++)
        for (j = 0; j < result.m.cols; ++j, result_i_j++)
            for (*result_i_j = k = 0; k < A->cols; k++)
                *result_i_j += *read_matrix_at(A, i, k) * *read_matrix_at(B, k, j);

    return result;
}

void cofactor_matrix_of(Matrix *A, int row, int col, Matrix *B)
{
    float *matB = B->mat;
    float *matA = A->mat;

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
    ResuM cof;
    
    if ((result.error = A->rows != A->cols))
        return result;

    result.f = 0;
    if (A->rows == 1)
        result.f = *A->mat;
    else {
        cof = init_matrix(A->rows - 1, A->cols - 1);
        if ((result.error = cof.error))
            return result;

        /* Tomar la primera fila y calcular el determinante de
         * los cofactores de cada entrada de la fila */
        for (int i = 0; i < A->rows; i++) {
            cofactor_matrix_of(A, 0, i, &cof.m);
            det_resu = determinant(&cof.m);

            if ((result.error = det_resu.error)) {
                free_matrix(&cof.m);
                return result;
            }

            if (i % 2 == 0)
                result.f += *(A->mat + i) * det_resu.f;
            else
                result.f -= *(A->mat + i) * det_resu.f;
        }

        free_matrix(&cof.m);
    }

    return result;
}

ResuM invert_matrix(Matrix *A)
{
    ResuF detA = determinant(A), detTemp;
    ResuM result = init_matrix(A->rows, A->cols),
          tempMat = init_matrix(A->rows - 1, A->cols - 1);
    int i, j;

    /* Return if det(A) == 0 */
    if ((result.error = detA.error || detA.f == 0))
        return result;

    if (result.error || tempMat.error) {
        if (!tempMat.error)
            free_matrix(&tempMat.m);
        if (!result.error)
            free_matrix(&result.m);

        result.error = 1;
        return result;
    }

    for (i = 0; i < A->rows; i++)
        for (j = 0; j < A->cols; j++) {
            cofactor_matrix_of(A, j, i, &tempMat.m);
            detTemp = determinant(&tempMat.m);

            if ((result.error = detTemp.error)) {
                free_matrix(&tempMat.m);
                return result;
            }

            *read_matrix_at(&result.m, i, j) =
                (((i + j) % 2 == 0) ? 1 : -1) * detTemp.f / detA.f;
        }

    free_matrix(&tempMat.m);
    return result;
}

ResuM transpose_matrix(Matrix *A)
{
    ResuM A_t = init_matrix(A->cols, A->rows);
    int i, j;

    if (!A_t.error)
        for (i = 0; i < A_t.m.rows; i++)
            for (j = 0; j < A_t.m.cols; j++)
                *read_matrix_at(&A_t.m, i, j)= *read_matrix_at(A, j, i);

    return A_t;
}

void print_matrix(Matrix *A)
{
    unsigned i, j;

    for (i = 0; i < A->rows; i++) {
        for (j = 0; j < A->cols; j++)
            printf("%10.3f ", *read_matrix_at(A, i, j));

        putchar('\n');
    }
}

void get_matrix(Matrix *A)
{
    float *mat;

    /* Scan rows * cols floats from stdin */
    for (mat = A->mat; mat - A->mat < A->rows * A->cols; mat++)
        scanf("%f", mat);
}
