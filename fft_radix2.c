#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_N 65536
#define PI 3.14159265358979323846

typedef struct {
    double re;
    double im;
} Complex;

static Complex add(Complex a, Complex b) {
    return (Complex){a.re + b.re, a.im + b.im};
}

static Complex sub(Complex a, Complex b) {
    return (Complex){a.re - b.re, a.im - b.im};
}

static Complex mul(Complex a, Complex b) {
    return (Complex){a.re*b.re - a.im*b.im,
                     a.re*b.im + a.im*b.re};
}

static int power_of_two(int n) {
    return n > 0 && (n & (n - 1)) == 0;
}

/* Forward, unnormalized, iterative radix-2 DIT FFT. */
static int fft(Complex *a, int n) {
    if (!a || !power_of_two(n) || n > MAX_N) return 0;
    if (n == 1) return 1;

    /* Incremental bit-reversed index: amortized O(n). */
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) {
            Complex temp = a[i];
            a[i] = a[j];
            a[j] = temp;
        }
    }

    for (int len = 2; ; len <<= 1) {
        int half = len / 2;
        for (int base = 0; base < n; base += len) {
            for (int j = 0; j < half; ++j) {
                double angle = -2.0 * PI * j / len;
                Complex w = {cos(angle), sin(angle)};
                Complex u = a[base + j];
                Complex v = mul(w, a[base + j + half]);
                a[base + j] = add(u, v);
                a[base + j + half] = sub(u, v);
            }
        }
        if (len == n) break;
    }
    return 1;
}

int main(void) {
    long requested;
    if (scanf("%ld", &requested) != 1) {
        fprintf(stderr, "ERROR: invalid N\n");
        return EXIT_FAILURE;
    }
    if (requested < 1 || requested > MAX_N ||
        !power_of_two((int)requested)) {
        fprintf(stderr, "ERROR: N must be a power of 2 (1..65536)\n");
        return EXIT_FAILURE;
    }
    int n = (int)requested;
    Complex *a = malloc((size_t)n * sizeof(*a));
    if (!a) {
        fprintf(stderr, "ERROR: allocation failed\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < n; ++i) {
        if (scanf("%lf %lf", &a[i].re, &a[i].im) != 2 ||
            !isfinite(a[i].re) || !isfinite(a[i].im) ||
            fabs(a[i].re) > 1e100 || fabs(a[i].im) > 1e100) {
            fprintf(stderr, "ERROR: invalid sample at index %d\n", i);
            free(a);
            return EXIT_FAILURE;
        }
    }
    if (!fft(a, n)) {
        fprintf(stderr, "ERROR: FFT failed\n");
        free(a);
        return EXIT_FAILURE;
    }
    puts("k real imag magnitude");
    for (int k = 0; k < n; ++k) {
        printf("%d %.9f %.9f %.9f\n", k,
               a[k].re, a[k].im, hypot(a[k].re, a[k].im));
    }
    free(a);
    return EXIT_SUCCESS;
}
