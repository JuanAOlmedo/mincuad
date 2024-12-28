#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "include/matrix.h"
#include "include/gc.h"

void exit_and_explain(void);

int main()
{
    Matrix A, A_transpose, coefmat, coefmat_inverse, Y, A_t_Y, X;
    int rows, cols, i, j;
    float x, y;
    atexit(GC_empty);

    printf("Ingrese grado de polinomio: ");
    scanf("%d", &cols);
    printf("Ingrese largo de la lista de datos: ");
    scanf("%d", &rows);

    if (rows < ++cols) {
        cols = rows;
        printf("INFO: Tomando un polinomio del grado del tamaño de la lista");
    }

    A = matrix_new(rows, cols);
    Y = matrix_new(rows, 1);

    for (i = 0; i < rows; i++) {
        scanf("%f %f", &x, &y);

        matrix_write(&Y, i, 0, y);
        for (j = 0; j < cols; j++)
            matrix_write(&A, i, j, (double) powf(x, (float) cols - j - 1));
    }

    A_transpose = matrix_transpose_of(A);
    coefmat = matrix_multiply(A_transpose, A);
    A_t_Y = matrix_multiply(A_transpose, Y);
    X = matrix_system_solve(coefmat, A_t_Y);

    puts("El polinomio que mejor aproxima la lista de datos es:");
    for (i = 0; i < X.rows; i++) {
        printf((matrix_read(X, i, 0) >= 0) ? "+ " : "- ");
        if (i < X.rows - 1)
            printf("%1.4f * x^%d ", fabs(matrix_read(X, i, 0)), X.rows - i - 1);
        else
            printf("%1.4f\n", fabs(matrix_read(X, i, 0)));
    }

    return 0;
}

void exit_and_explain(void)
{
    fprintf(stderr, "Error de memoria\n");
    exit(EXIT_FAILURE);
}
