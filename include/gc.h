#ifndef _GC_H
#define _GC_H
#include "matrix.h"
#include <stdlib.h>

/* Añadir variable a lista */
void GC_push(void *var, void free_fun(void *));

/* Sacar variable de lista */
void GC_remove(void *var);

void GC_empty(void);

#endif
