#ifndef _MATRIX_H
#define _MATRIX_H

typedef struct {
   unsigned rows, cols;
   double *mat;
} Matrix;

struct LU {
    Matrix L, U;
    unsigned *p, n;
};

/* Devuelve un Matrix que contiene la matriz de tamaño rows x cols. */
Matrix matrix_new(unsigned rows, unsigned cols);

/* Devuelve 0 cuando la matriz no es válida */
int matrix_is_valid(Matrix A);

/* Devuelve una copia de A */
Matrix matrix_copy(Matrix A);

/* Devuelve A[row, col] */
double matrix_read(Matrix A, unsigned row, unsigned col);

/* Escribe x en A[row, col] */
void matrix_write(Matrix *A, unsigned row, unsigned col, float x);

/* Devuelve un arreglo con la columna n-ésima de A */
double *matrix_col(Matrix A, unsigned col);

/* Devuelve un arreglo con la fila n-ésima de A */
double *matrix_row(Matrix A, unsigned row);

/* Libera la memoria asignada a la matriz A */
void matrix_free(Matrix *A);

/* Devuelve A * x */
Matrix matrix_smultiply(Matrix A, double x);

/* Devuelve A * B. Devuelve una matriz inválida si las matrices son
 * no conformantes */
Matrix matrix_multiply(Matrix A, Matrix B);

/* Devuelve A + B. Devuelve una matriz inválida si las matrices son
 * no conformantes */
Matrix matrix_add(Matrix A, Matrix B);

/* Imprime A a la salida estándar */
void matrix_print(Matrix A);

/* Lee A de la entrada estándar */
void matrix_get(Matrix *A);

/* Devuelve det(A) si la matriz es cuadrada. Si no, devuelve 0 */
double matrix_det(Matrix A);

/* Devuelve la transpuesta de A */
Matrix matrix_transpose_of(Matrix A);

/* Transpone A */
void matrix_transpose(Matrix *A);

/* Devuelve la inversa de A si existe, un matriz inválida si no */
Matrix matrix_inv(Matrix A);

/* Devuelve la descomposición LU de A */
struct LU matrix_lu(Matrix A);

/* Devuelve la solución del sistema A * x = b */
Matrix matrix_system_solve(Matrix A, Matrix b);

#endif
