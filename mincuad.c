#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gc.h"

static ResuM A, A_transpose, coefmat, coefmat_inverse, Y, A_t_Y, X;

/* Devuelve n^m */
int expon(int n, int m);

void exit_and_explain(void);

int main()
{/*
    int rows, cols;
    char arr[100];
    scanf("%d", &rows);
    scanf("%d", &cols);
    A = init_matrix(rows, cols);
    scanf("%d", &rows);
    scanf("%d", &cols);
    Y = init_matrix(rows, cols);

    get_matrix(&A.m);
    get_matrix(&Y.m);

    printf("> ");
    scanf("%s", arr);
    while (strcmp(arr, "End") != 0) {
        if (strcmp(arr, "Transpose") == 0) {
            A_transpose = transpose_matrix(&A.m);
            print_matrix(&A_transpose.m);
        } else if (strcmp(arr, "Multiply") == 0) {
            A_transpose = multiply_matrix(&A.m, &Y.m);
            print_matrix(&A_transpose.m);
        } else if (strcmp(arr, "Determinant") == 0) {
            ResuF det = determinant(&A.m);
            printf("%f\n", det.f);
        } else if (strcmp(arr, "Invert") == 0) {
            A_transpose = invert_matrix(&A.m);
            if (A_transpose.error)
                puts("Error");
            print_matrix(&A_transpose.m);
        } else if (strcmp(arr, "Print") == 0) {
            print_matrix(&A.m);
            print_matrix(&Y.m);
        }

        printf("> ");
        scanf("%s", arr);
    } 
    */
    int rows, cols, i, j, x, y;
    atexit(free_all);

    printf("Ingrese grado de polinomio: ");
    scanf("%d", &cols);
    cols++;
    printf("Ingrese largo de la lista de datos: ");
    scanf("%d", &rows);

    A = init_matrix(rows, cols);
    Y = init_matrix(rows, 1);

    if (A.error || Y.error)
        exit_and_explain();

    for (i = 0; i < rows; i++) {
        scanf("%d %d", &x, &y);

        *read_matrix_at(&Y.m, i, 0) = y;
        for (j = 0; j < cols; j++)
            *read_matrix_at(&A.m, i, j) = expon(x, cols - j - 1);
    }

    A_transpose = transpose_matrix(&A.m);
    if (A_transpose.error)
        exit_and_explain();

    coefmat = multiply_matrix(&A_transpose.m, &A.m);
    A_t_Y = multiply_matrix(&A_transpose.m, &Y.m);
    free_matrix(&A.m);
    free_matrix(&A_transpose.m);
    free_matrix(&Y.m); 

    if (coefmat.error || A_t_Y.error)
        exit_and_explain();

    coefmat_inverse = invert_matrix(&coefmat.m);
    free_matrix(&coefmat.m);

    if (coefmat_inverse.error)
        exit_and_explain();

    X = multiply_matrix(&coefmat_inverse.m, &A_t_Y.m);
    free_matrix(&coefmat_inverse.m);
    free_matrix(&A_t_Y.m);

    if (X.error)
        exit_and_explain();

    print_matrix(&X.m);

    return 0;
}

int expon(int n, int m)
{
    if (m == 0)
        return 1;

    int resu = 1;
    for (int i = 0; i < m; i++)
        resu *= n;

    return resu;
}

void exit_and_explain(void)
{
    fprintf(stderr, "Error de memoria\n");
    exit(EXIT_FAILURE);
}
