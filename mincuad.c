#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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

    A = init_matrix(rows, cols);
    Y = init_matrix(rows, 1);

    if (A.mat == NULL || Y.mat == NULL)
        exit_and_explain();

    for (i = 0; i < rows; i++) {
        scanf("%f %f", &x, &y);

        *read_matrix_at(&Y, i, 0) = y;
        for (j = 0; j < cols; j++)
            *read_matrix_at(&A, i, j) =(double) powf(x, (float) cols - j - 1);
    }

    A_transpose = transpose_matrix_of(A);
    if (A_transpose.mat == NULL)
        exit_and_explain();

    coefmat = multiply_matrices(&A_transpose, &A);
    A_t_Y = multiply_matrices(&A_transpose, &Y);

    if (coefmat.mat == NULL || A_t_Y.mat == NULL)
        exit_and_explain();

    X = solve(coefmat, A_t_Y);

    if (X.mat == NULL)
        exit_and_explain();

    puts("El polinomio que mejor aproxima la lista de datos es:");
    for (i = 0; i < X.rows; i++) {
        printf((*read_matrix_at(&X, i, 0) >= 0) ? "+ " : "- ");
        if (i < X.rows - 1)
            printf("%1.4f * x^%d ", fabs(*read_matrix_at(&X, i, 0)), X.rows - i - 1);
        else
            printf("%1.4f\n", fabs(*read_matrix_at(&X, i, 0)));
    }

    return 0;
}

void exit_and_explain(void)
{
    fprintf(stderr, "Error de memoria\n");
    exit(EXIT_FAILURE);
}
