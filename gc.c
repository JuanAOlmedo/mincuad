#include <stdlib.h>
#include "include/matrix.h"

static struct matrix_array {
    struct matrix_array *nxt;
    double *mat;
} *matrices = NULL;

void remove_from_matrix_list(double *mat)
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

void push_to_matrix_list(Matrix *A)
{
    struct matrix_array *it, *el = malloc(sizeof(struct matrix_array));    

    if (el == NULL)
        return;
    
    if (matrices == NULL)
        matrices = el;
    else {
        for (it = matrices; it->nxt != NULL; it = it->nxt)
            ;
        it->nxt = el;
    }

    el->mat = A->mat;
    el->nxt = NULL;
}

void free_all(void)
{
    double *aux;

    while (matrices != NULL) {
        aux = matrices->mat;
        remove_from_matrix_list(aux);
        free(aux);
    }
}
