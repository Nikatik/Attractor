#include <math.h>
#include <quadmath.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define EPS (__float128) pow (10, -34)
#define pi  3.1415926535897932384626433

__float128 astep (
    __float128 T, __float128* x, __float128* y, __float128* z,
    long long unsigned* i, long long unsigned* j, unsigned int p,
    unsigned int s, __float128** k, __float128** cab, __float128 tol,
    __float128 f (__float128, __float128, __float128, __float128, __float128),
    __float128 g (__float128, __float128, __float128, __float128, __float128),
    __float128 q (__float128, __float128, __float128, __float128, __float128),
    __float128 R, bool first);
