#include "matrix.h"

/* Acumula todas las matrices definidas en una lista
 * para ser liberadas al final del programa */

/* Libera todas las matrices declaradas */
void free_all(void);

/* Añadir matriz a lista */
int push_to_matrix_list(Matrix *A);

/* Sacar matriz de lista */
void remove_from_matrix_list(double *mat);
