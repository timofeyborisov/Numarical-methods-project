#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

enum {FUNC = 0, POLY = 1, SPLN = 2};

int num = 101;
// number of dots for graph
int coeff = 1;
// coefficient (accuracy) for max deviation calculation

typedef struct {
    double a, b, c, d;
} Spline;

int
factorial(int n)
{
    int ans = 1;
    for (int i = 1; i <= n; ++i) {
        ans *= i;
    }
    return ans;
}

double
fifth_derivative(double x)
{
    // d5 / dx5 (x^2 - 1)^5
    return 480 * x * (63 * pow(x, 4) - 70 * pow(x, 2) + 15);
}

double
f1(double x)
{
    // first function
    // first data set

    return fifth_derivative(x) / (pow(2, 5) * factorial(5));
}

double
f2(double x)
{
    // second data set
    return sin(4 * x) * exp(2 * x);
}

double 
oriented_polynomial(double x, double x0, int i, int n, double h)
{
    // x - the value of the oriented polynomial at point x
    // x0 - left boundary of a set (segment)
    // i - number of xi
    // n - number of interpolation nodes
    // h - grid step

    double ans = 1;
    double xi = x0 + i * h;
    for (int j = 0; j < n; ++j) {
        if (j == i) {
            continue;
        }
        double xj = x0 + j * h;
        ans *= (x - xj) / (xi - xj);
    }
    return ans;
}

double
Lagrange_polynomial(double (*f)(double), double x, double x0, double xn, int n)
{
    // f - the function for which the Lagrange polynomial is calculated
    // x - the value of the Lagrange polynomial at point x
    // x0 - left boundary of a set (segment)
    // xn - right boundary of a set 
    // n - number of interpolation nodes

    double ans = 0;
    double h = (xn - x0) / (n - 1);
    for (int i = 0; i < n; ++i) {
        double xi = x0 + i * h;
        ans += f(xi) * oriented_polynomial(x, x0, i, n, h);
    }
    return ans;
}

double
spline_val(Spline coeff, double x, double xi)
{
    return (coeff.a + coeff.b * (x - xi) + coeff.c / 
    2 * pow(x - xi, 2) + coeff.d / 6 * pow(x - xi, 3));
}

double
fi(double (*f)(double), double xi, double h)
{
    return 6 * (f(xi - h) - 2 * f(xi) + f(xi + h)) / pow(h, 2);
}

void
run_method_c(double (*f)(double), Spline *coeff_arr, double x0, int n, double h)
{
    // n - number of n - number of segments
    // h - step

    double *alpha = malloc(n * sizeof *alpha);
    double *beta = malloc(n * sizeof *beta);
    double A = 1, C = 4, B = 1;
    // A, B, C - const

    coeff_arr[0].c = 0;
    coeff_arr[n].c = 0;
    
    alpha[1] = 0;
    beta[1] = 0;
    
    for (int i = 1; i <= n - 1; ++i) {
        double xi = x0 + i * h;
        // i = 1 .. n - 1
        double denom = A * alpha[i] + C;
        double F = fi(f, xi, h);

        alpha[i + 1] = -B / denom;
        beta[i + 1] = (F - A * beta[i]) / denom;
    }

    for (int i = n - 1; i >= 0; --i) {
        coeff_arr[i].c = alpha[i + 1] * coeff_arr[i + 1].c + beta[i + 1];
    }

    free(alpha);
    free(beta);
}

void
get_coeff_abd(double (*f)(double), Spline *coeff_arr, double x0, int n, double h)
{
    // f - the function for which the system solution is calculated
    // x0 - left boundary of a set (segment)
    // xn - right boundary of a set 
    // n - number of segments

    for (int i = 1; i <= n; ++i) {
        double xi = x0 + i * h;
        coeff_arr[i].a = f(xi);
        coeff_arr[i].d = (coeff_arr[i].c - coeff_arr[i - 1].c) / h;
        coeff_arr[i].b = h * coeff_arr[i].c / 2 - pow(h, 2) * coeff_arr[i].d / 6 + 
        (f(xi) - f(xi - h)) / h;
    }
}

double
cubic_spline(double (*f)(double), double x, double x0, double xn, int n)
{
    // f - the function for which the cubic spline S(x) is calculated
    // x - the value of the cubic spline Si(x) at point x
    // x0 - left boundary of a set (segment)
    // xn - right boundary of a set 
    // n - number of splines (segments)

    Spline *coeff_arr = malloc((n + 1) * sizeof *coeff_arr);
    // coeff_arr[i] - i-coefficients of Si

    double h = (xn - x0) / n;
    // step
    run_method_c(f, coeff_arr, x0, n, h);
    // calculating c
    get_coeff_abd(f, coeff_arr, x0, n, h);
    // calculating a, b and d using c
    
    double ans = 0;
    for (int i = 1; i <= n; ++i) {
        double l = x0 + (i - 1) * h;
        // xi - 1
        double r = x0 + i * h;
        // xi
        if (l <= x && x <= r) {
            ans = spline_val(coeff_arr[i], x, r);
            break;
        }
    }

    free(coeff_arr);
    return ans;
}

double
max_deviation(double (*f)(double), int type, double x0, double xn, int points, int n)
{
    double ans = 0;
    double h = (xn - x0) / (points - 1);
    double y = 0;
    for (int i = 0; i < points; ++i) {
        double x = x0 + i * h;
        switch (type) {
            case POLY: {
                y = Lagrange_polynomial(f, x, x0, xn, n);
                break;
            }
            case SPLN: {
                y = cubic_spline(f, x, x0, xn, n - 1);
                break;
            }
        }
        double diff = fabs(f(x) - y);
        if (diff > ans) {
            ans = diff;
        }
    }
    return ans;
}

void
make_data_file(char *filename, double (*f)(double), int type, double x0, double xn, int n)
{
    // filename - name of created file
    // f - function
    // type:
    // FUNC - function
    // POLY - calculating Largrange polynomial for function
    // SPLN - spline interpolation for function 

    int fd = open(filename, O_CREAT | O_TRUNC | O_WRONLY, 0666);
    char buf[256];
    double h = (xn - x0) / (num - 1);
    double y = 0;
    for (int i = 0; i < num; ++i) {
        double x = x0 + i * h;
        switch (type) {
            case FUNC: {
                y = f(x);
                break;
            }
            case POLY: {
                y = Lagrange_polynomial(f, x, x0, xn, n);
                break;
            }
            case SPLN: {
                y = cubic_spline(f, x, x0, xn, n - 1);
                break;
            }
            default:
        }
        snprintf(buf, sizeof buf, "%lf %lf\n", x, y);
        write(fd, buf, strlen(buf));
    }
    close(fd);
}

int
main(void)
{
    printf("Enter number of dots for graphs: ");
    int res = scanf("%d", &num);
    if (res == -1 || num <= 2) {
        printf("Wrong value\n");
        return 0;
    }

    make_data_file("f1.dat", f1, FUNC, 0, 2, 0);
    make_data_file("f2.dat", f2, FUNC, 0, 2, 0);
    // for f1, f2

    make_data_file("f1_poly_3.dat", f1, POLY, 0, 2, 3);
    make_data_file("f1_poly_5.dat", f1, POLY, 0, 2, 5);
    make_data_file("f1_poly_9.dat", f1, POLY, 0, 2, 9);
    make_data_file("f1_poly_17.dat", f1, POLY, 0, 2, 17);
    // for Lagrange polynomials (for f1)

    make_data_file("f2_poly_3.dat", f2, POLY, 0, 2, 3);
    make_data_file("f2_poly_5.dat", f2, POLY, 0, 2, 5);
    make_data_file("f2_poly_9.dat", f2, POLY, 0, 2, 9);
    make_data_file("f2_poly_17.dat", f2, POLY, 0, 2, 17);
    // for Lagrange polynomials (for f2)
    
    make_data_file("f2_spln_3.dat", f2, SPLN, 0, 2, 3);
    make_data_file("f2_spln_5.dat", f2, SPLN, 0, 2, 5);
    make_data_file("f2_spln_9.dat", f2, SPLN, 0, 2, 9);
    make_data_file("f2_spln_17.dat", f2, SPLN, 0, 2, 17);
    // for cubic splines

    printf("---\n");
    printf("All data files were successfully created\n");
    
    double max_dev3 = max_deviation(f1, POLY, 0, 2, 1001, 3);
    double max_dev5 = max_deviation(f1, POLY, 0, 2, 1001, 5);
    double max_dev9 = max_deviation(f1, POLY, 0, 2, 1001, 9);
    double max_dev17 = max_deviation(f1, POLY, 0, 2, 1001, 17);

    printf("---\n");
    printf("Max deviation for f1 and Lagrange polynomial (n = 3): %lf\n", max_dev3);
    printf("Max deviation for f1 and Lagrange polynomial (n = 5): %lf\n", max_dev5);
    printf("Max deviation for f1 and Lagrange polynomial (n = 9): %lf\n", max_dev9);
    printf("Max deviation for f1 and Lagrange polynomial (n = 17): %lf\n", max_dev17);

    max_dev3 = max_deviation(f2, POLY, 0, 2, 1001, 3);
    max_dev5 = max_deviation(f2, POLY, 0, 2, 1001, 5);
    max_dev9 = max_deviation(f2, POLY, 0, 2, 1001, 9);
    max_dev17 = max_deviation(f2, POLY, 0, 2, 1001, 17);

    printf("---\n");
    printf("Max deviation for f2 and cubic spline (n = 3): %lf\n", max_dev3);
    printf("Max deviation for f2 and cubic spline (n = 5): %lf\n", max_dev5);
    printf("Max deviation for f2 and cubic spline (n = 9): %lf\n", max_dev9);
    printf("Max deviation for f2 and cubic spline (n = 17): %lf\n", max_dev17);

    printf("---\n");
    return 0;
}