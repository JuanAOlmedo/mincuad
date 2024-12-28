#ifndef _GC_H
#define _GC_H
#include <stdlib.h>

/* Acumula todas las variables definidas en una lista
 * para ser liberadas al final del programa */

/* Añadir variable a lista */
void GC_push(void *var, void free_fun(void *));

/* Sacar variable de lista */
void GC_remove(void *var);

/* Vaciar lista */
void GC_empty(void);

#endif
