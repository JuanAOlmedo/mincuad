#ifndef _GC_H
#define _GC_H
#include "matrix.h"
#include <stdlib.h>

/* Acumula todas las matrices definidas en una lista
 * para ser liberadas al final del programa */

struct matrix_list {
    struct matrix_list *nxt;
    double *mat;
};

/* Devuelve la lista en donde se guardan las matrices */
struct matrix_list *get_matrix_list(void);

/* Añadir matriz a lista */
void push_to_matrix_list(Matrix *A);

/* Sacar matriz de lista */
void remove_from_matrix_list(double *mat);

#endif
