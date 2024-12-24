#ifndef _MATRIX_H
#define _MATRIX_H

typedef struct {
   unsigned rows, cols;
   double *mat;
} Matrix;

typedef struct {
    double f;
    short error;
} ResuF;

struct LU {
    Matrix L, U;
    unsigned *p, n;
};

/* Devuelve un Matrix que contiene la matriz de tamaño rows x cols.
 * Si no hay memoria suficiente, devuelve en estado de error. */
Matrix init_matrix(unsigned rows, unsigned cols);

/* Devuelve una copia de A */
Matrix copy_matrix(Matrix A);

/* Devuelve M[row, col] */
double *read_matrix_at(Matrix *M, unsigned row, unsigned col);

/* Equivalente a read_matrix_at(M, 0, col) */
double *col(Matrix *M, unsigned col);

/* Equivalente a read_matrix_at(M, row, 0) */
double *row(Matrix *M, unsigned row);

/* Libera el puntero a la matriz */
void free_matrix(Matrix *M);

/* Devuelve A * x. Devuelve en estado de error si no hay memoria
 * disponible */
Matrix multiply_matrix(Matrix *A, double x);

/* Devuelve A * B. Devuelve en estado de error si las matrices son
 * no conformantes o si no hay memoria disponible */
Matrix multiply_matrices(Matrix *A, Matrix *B);

/* Devuelve A + B. Devuelve en estado de error si las matrices son
 * no conformantes o si no hay memoria disponible */
Matrix add_matrices(Matrix *A, Matrix *B);

/* Imprime A */
void print_matrix(Matrix *A);

/* Lee A de la entrada estándar */
void get_matrix(Matrix *A);

/* Devuelve det(A), error si la matrix no es cuadrada */
ResuF determinant(Matrix *A);

/* Devuelve la transpuesta de A */
Matrix transpose_matrix_of(Matrix A);

/* Escribe en A su transpuesta */
void transpose_matrix(Matrix *A);

/* Devuelve la inversa de A si existe, error si no o si no hay
 * memoria suficiente */
Matrix invert_matrix(Matrix *A);

/* Devuelve la descomposición LU de A */
struct LU lu(Matrix A);
#endif
