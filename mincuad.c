#include <stdio.h>
#include <math.h>
#include "include/matrix.h"

int approx(double x, double y)
{
	return fabs(x - y) < 1e-8;
}

int main()
{
    Matrix A, t, y, x;
    int rows = -1, cols = -1;

	while (20 <= cols || cols < 0) {
		printf("Ingrese grado de polinomio (0-20): ");
		scanf("%d", &cols);
	}

	// El número de columnas de la matriz es 1 + el grado del polinomio
	cols++;

	while (rows < cols) {
		printf("Ingrese largo de la lista de datos (%d-): ", cols);
		scanf("%d", &rows);
	}

    t = matrix_new(rows, 1);
    y = matrix_new(rows, 1);

	puts("Ingrese x:");
    matrix_get(&t);

    A = matrix_vander(t, cols);

	puts("Ingrese y:");
    matrix_get(&y);

    x = matrix_ls_solve(A, y);

    printf("El polinomio que mejor aproxima la lista de datos es:\np(x) = ");
    for (unsigned i = 0; i < x.rows; i++) {
		double x_i = matrix_read(x, i, 0);
		// Mostrar el término cuando el coeficiente no es 0
		if (approx(fabs(x_i), 1)) {
			// Si x_i == 1 o x_i == -1, no mostrar el coeficiente exacto
			printf((x_i >= 0) ? "+ " : "- ");
			printf("x^%u ", x.rows - i - 1);
		} else if (!approx(x_i, 0)) {
			printf((x_i >= 0) ? "+ " : "- ");
			if (i < x.rows - 1)
				printf("%1.7fx^%u ", fabs(x_i), x.rows - i - 1);
			else
				printf("%1.7f", fabs(x_i)); // Término constante
		}
    }
	putchar('\n');
	
    matrix_free(&x);
    matrix_free(&A);
    matrix_free(&t);
    matrix_free(&y);

    return 0;
}
