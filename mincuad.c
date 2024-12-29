#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "include/matrix.h"
#include "include/gc.h"

void exit_and_explain(void);

int main()
{
    Matrix A, A_t, A_t_A, T, Y, A_t_Y, X;
    int rows, cols;
    atexit(GC_empty);

    printf("Ingrese grado de polinomio: ");
    scanf("%d", &cols);
    printf("Ingrese largo de la lista de datos: ");
    scanf("%d", &rows);

    if (rows < ++cols) {
        cols = rows;
        printf("INFO: Tomando un polinomio del grado del tamaño de la lista");
    }

    T = matrix_new(rows, 1);
    Y = matrix_new(rows, 1);

    matrix_get(&T);
    matrix_get(&Y);

    A = matrix_vander(T, cols);

    A_t = matrix_transpose_of(A);
    A_t_A = matrix_multiply(A_t, A);
    A_t_Y = matrix_multiply(A_t, Y);
    X = matrix_system_solve(A_t_A, A_t_Y);

    puts("El polinomio que mejor aproxima la lista de datos es:");
    for (unsigned i = 0; i < X.rows; i++) {
        printf((matrix_read(X, i, 0) >= 0) ? "+ " : "- ");
        if (i < X.rows - 1)
            printf("%1.4f * x^%d ", fabs(matrix_read(X, i, 0)), X.rows - i - 1);
        else
            printf("%1.4f\n", fabs(matrix_read(X, i, 0)));
    }
    matrix_free(&A);
    matrix_free(&T);
    matrix_free(&Y);
    matrix_free(&X);
    matrix_free(&A_t);
    matrix_free(&A_t_Y);
    matrix_free(&A_t_A);

    return 0;
}

void exit_and_explain(void)
{
    fprintf(stderr, "Error de memoria\n");
    exit(EXIT_FAILURE);
}
