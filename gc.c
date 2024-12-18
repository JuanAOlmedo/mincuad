#include "include/gc.h"

/* Acumula todas las variables definidas en una lista
 * para ser liberadas al final del programa */
struct variable_list {
    struct variable_list *nxt;
    void *var;
    void (*free_fun)(void *);
};

static struct variable_list *variables = NULL;

void GC_remove(void *var)
{
    struct variable_list *aux, *it = variables;

    if (variables->var != var) {
        while (it->nxt != NULL && it->nxt->var != var)
            it = it->nxt;

        aux = it->nxt;
        it->nxt = aux->nxt;
    } else {
        aux = variables;
        variables = variables->nxt;
    }

    aux->free_fun(aux->var);
    free(aux);
}

void GC_push(void *var, void (*free_fun)(void *))
{
    struct variable_list *it = variables,
        *el = malloc(sizeof(struct variable_list));

    if (el == NULL)
        return;

    if (variables != NULL) {
        while (it->nxt != NULL)
            it = it->nxt;
        it->nxt = el;
    } else
        variables = el;

    *el = (struct variable_list) {NULL, var, free_fun};
}

void GC_empty(void)
{
    struct variable_list *aux;

    while (variables != NULL) {
        aux = variables;
        variables = variables->nxt;
        aux->free_fun(aux->var);
        free(aux);
    }
}
