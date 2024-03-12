typedef struct {
   unsigned rows, cols;
   float *mat;
} Matrix;

typedef struct {
    float f;
    short error;
} ResuF;

/* Devuelve un Matrix que contiene la matriz de tamaño rows x cols.
 * Si no hay memoria suficiente, devuelve en estado de error. */
Matrix init_matrix(unsigned rows, unsigned cols);

/* Devuelve M[row, col] */
float *read_matrix_at(Matrix *M, unsigned row, unsigned col);

/* Equivalente a read_matrix_at(M, 0, col) */
float *col(Matrix *M, unsigned col);

/* Equivalente a read_matrix_at(M, row, 0) */
float *row(Matrix *M, unsigned row);

/* Libera el puntero a la matriz */
void free_matrix(Matrix *M);

/* Devuelve A * B. Devuelve en estado de error si las matrices son
 * no conformantes o si no hay memoria disponible */
Matrix multiply_matrix(Matrix *A, Matrix *B);

/* Imprime A */
void print_matrix(Matrix *A);

/* Lee A de la entrada estándar */
void get_matrix(Matrix *A);

/* Devuelve det(A), error si la matrix no es cuadrada */
ResuF determinant(Matrix *A);

/* Devuelve la transpuesta de A, error si no hay memoria suficiente */
Matrix transpose_matrix(Matrix *A);

/* Devuelve la inversa de A si existe, error si no o si no hay
 * memoria suficiente */
Matrix invert_matrix(Matrix *A);
