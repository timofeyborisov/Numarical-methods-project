#include <stdio.h>
#include <math.h>

double 
f(double x)
{
    return (pow(x * x - 1, 5));
}

double 
der(double (*f)(double), double x, double h)
{
    return (f(x - 3 * h) - 5 * f(x - 2 * h) + 10 * f(x - h) - 10 * f(x + h) + 5 * f(x + 2 * h) - f(x + 3 * h));
}

double 
f1(double x, double h)
{
    return x / 3840;
}

int main() {
    double x;
    double h = 0.01;
    scanf("%lf", &x);

    printf("%lf\n", der(f, x, h));
    return 0;
}