#include "include/gc.h"
#include "include/matrix.h"

static struct matrix_list *matrices = NULL;

struct matrix_list *get_matrix_list(void)
{
    return matrices;
}

void remove_from_matrix_list(double *mat)
{
    struct matrix_list *aux, *it = matrices;

    if (matrices->mat != mat) {
        for (; it->nxt != NULL && it->nxt->mat != mat; it = it->nxt)
            ;

        aux = it->nxt;
        it->nxt = aux->nxt;
    } else {
        aux = matrices;
        matrices = matrices->nxt;
    }

    free(aux);
}

void push_to_matrix_list(Matrix *A)
{
    struct matrix_list *it, *el = malloc(sizeof(struct matrix_list));    

    if (el == NULL)
        return;
    
    if (matrices == NULL)
        matrices = el;
    else {
        for (it = matrices; it->nxt != NULL; it = it->nxt)
            ;
        it->nxt = el;
    }

    *el = (struct matrix_list) {NULL, A->mat};
}

