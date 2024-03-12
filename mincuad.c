#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "include/gc.h"

static ResuM A, A_transpose, coefmat, coefmat_inverse, Y, A_t_Y, X;

void exit_and_explain(void);

int main()
{
    int rows, cols, i, j;
    float x, y;
    atexit(free_all);

    printf("Ingrese grado de polinomio: ");
    scanf("%d", &cols);
    printf("Ingrese largo de la lista de datos: ");
    scanf("%d", &rows);

    if (rows > ++cols) {
        cols = rows;
        printf("INFO: Tomando un polinomio del grado del tamaño de la lista");
    }

    A = init_matrix(rows, cols);
    Y = init_matrix(rows, 1);

    if (A.error || Y.error)
        exit_and_explain();

    for (i = 0; i < rows; i++) {
        scanf("%f %f", &x, &y);

        *read_matrix_at(&Y.m, i, 0) = y;
        for (j = 0; j < cols; j++)
            *read_matrix_at(&A.m, i, j) = powf(x, (float) cols - j - 1);
    }

    A_transpose = transpose_matrix(&A.m);
    if (A_transpose.error)
        exit_and_explain();

    coefmat = multiply_matrix(&A_transpose.m, &A.m);
    A_t_Y = multiply_matrix(&A_transpose.m, &Y.m);

    if (coefmat.error || A_t_Y.error)
        exit_and_explain();

    coefmat_inverse = invert_matrix(&coefmat.m);

    if (coefmat_inverse.error)
        exit_and_explain();

    X = multiply_matrix(&coefmat_inverse.m, &A_t_Y.m);

    if (X.error)
        exit_and_explain();

    print_matrix(&X.m);

    return 0;
}

void exit_and_explain(void)
{
    fprintf(stderr, "Error de memoria\n");
    exit(EXIT_FAILURE);
}
