#include <stdio.h>
#include <math.h>
#include "include/matrix.h"

int main()
{
    Matrix A, t, y, x;
    int rows, cols;

    printf("Ingrese grado de polinomio: ");
    scanf("%d", &cols);
    printf("Ingrese largo de la lista de datos: ");
    scanf("%d", &rows);

    if (rows < ++cols) {
        cols = rows;
        printf("INFO: Tomando un polinomio del grado del tamaño de la lista");
    }

    t = matrix_new(rows, 1);
    y = matrix_new(rows, 1);

    matrix_get(&t);
    matrix_get(&y);

    A = matrix_vander(t, cols);
    matrix_free(&t);

    x = matrix_ls_solve(A, y);
    matrix_free(&A);
    matrix_free(&y);

    puts("El polinomio que mejor aproxima la lista de datos es:");
    for (unsigned i = 0; i < x.rows; i++) {
        printf((matrix_read(x, i, 0) >= 0) ? "+ " : "- ");
        if (i < x.rows - 1)
            printf("%1.6f * x^%d ", fabs(matrix_read(x, i, 0)), x.rows - i - 1);
        else
            printf("%1.6f\n", fabs(matrix_read(x, i, 0)));
    }
    matrix_free(&x);

    return 0;
}
