#include "lib.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

static __float128 f (__float128 t, __float128 x, __float128 y, __float128 z,
                     __float128 R)        //   y` = f(t,x,y,z)
{
    return 10.Q * (y - x);
}

static __float128 g (__float128 t, __float128 x, __float128 y, __float128 z,
                     __float128 R)        //   x` = g(t,x,y,z)
{
    return R * x - y - x * z;
}

static __float128 q (__float128 t, __float128 x, __float128 y, __float128 z,
                     __float128 R)        //   z` = q(t,x,y,z)
{
    return -8.Q / 3 * z + x * y;
}

#pragma GCC diagnostic pop

int mainq (unsigned int, __float128);
int open (FILE**, char*, unsigned int*, unsigned int*);

__float128 rd (FILE*);

void initialization (long long unsigned, long long unsigned, unsigned int,
                     __float128***, __float128***);
void reading (FILE*, __float128**, unsigned int, long long unsigned,
              long long unsigned);
void printing (__float128**, unsigned int, long long unsigned,
               long long unsigned);
void cleaning (long long unsigned, unsigned int, FILE*, __float128**,
               __float128**);

int main (int argc, char* argv[])
{
    if (argc > 1)
    {
        unsigned int n;
        long double  R = 1.Q;
        sscanf (argv[1], "%u", &n);
        if (argc > 2) { sscanf (argv[2], "%Lf", &R); }
        return mainq (n, (__float128) R);
    }
    return mainq (12, (__float128) 1);
}

int mainq (unsigned int toller, __float128 R)
{
    __float128 tol = (__float128) pow (10, -(int) toller);

    __float128 diff = 0.00000005Q;

    __float128 x = diff;
    __float128 y = diff;
    __float128 z = 0;
    FILE*      inpf;

    __float128 T = (__float128) 1000 * powq(R/4 , -5.Q/3);
    __float128 err;

    unsigned int s;
    unsigned int p;

    __float128** k;
    __float128** cab;

    clock_t timer;

    long long unsigned i = 0;
    long long unsigned j = 0;

    /////////////////////////////////////////////////////////////////////////////////////////

    switch (open (&inpf, "koef (8).txt", &s, &p))
    {
        case -1: return -1;
        case -2: return -2;
        default: printf ("File opened.\n");
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    initialization (i, j, s, &cab, &k);
    printf ("All initialized.\n");

    /////////////////////////////////////////////////////////////////////////////////////////

    reading (inpf, cab, s, i, j);
    printf ("RK matrix readed.\n\n");

    // printing (cab, s, i, j);

    /////////////////////////////////////////////////////////////////////////////////////////

    // The shoting method
    timer = clock();

    err = astep (T, &x, &y, &z, &i, &j, p, s, k, cab, tol, f, g, q, R, true);
    printf ("The Runge-Kutta output:\n%.10Le    %.10Le    %.10Le  |  %.2Ld    "
            "%.2Ld  | "
            " %.0Le    %.7Le  |  %.3Le    %.3Le\n\n",
            (long double) x, (long double) y, (long double) z, i, j,
            (long double) tol, (long double) err, (long double) T,
            (long double) R);


    // x = -diff;
    // y = -diff;
    // z = 0;
    // i = j = 0;

    // err = astep (T, &x, &y, &z, &i, &j, p, s, k, cab, tol, f, g, q, R, false);
    // printf ("The Runge-Kutta output:\n%.10Le    %.10Le    %.10Le  |  %.2Ld    "
    //         "%.2Ld  | "
    //         " %.0Le    %.7Le  |  %.3Le    %.3Le\n\n",
    //         (long double) x, (long double) y, (long double) z, i, j,
    //         (long double) tol, (long double) err, (long double) T,
    //         (long double) R);

    timer -= clock();
    printf ("%.5Lf seconds\n", ((long double) -timer) / CLOCKS_PER_SEC);

    /////////////////////////////////////////////////////////////////////////////////////////

    cleaning (i, s, inpf, cab, k);
    printf ("All cleaned.\n");

    /////////////////////////////////////////////////////////////////////////////////////////

    return 0;
}

// File checker
int open (FILE** inpf, char* path, unsigned int* s, unsigned int* p)
{

    *inpf = fopen (path, "r");
    if (*inpf == NULL)
    {
        printf ("File doen`t exist\n");
        return -1;
    }
    if (!fscanf (*inpf, "%ud", p) || !fscanf (*inpf, "%ud", s))
    {
        printf ("File isn`t correct\n");
        fclose (*inpf);
        return -2;
    }
    return 0;
}

// Constructor
void initialization (long long unsigned i, long long unsigned j, unsigned int s,
                     __float128*** cab, __float128*** k)
{
    *cab = (__float128**) malloc ((s + 3) * sizeof (__float128*));
    for (i = 0; i < s + 3; i++)
    {
        (*cab)[i] = (__float128*) malloc (s * sizeof (__float128));
        for (j = 0; j < s; j++) (*cab)[i][j] = 0;
    }

    *k = (__float128**) malloc (3 * sizeof (__float128*));
    for (i = 0; i < 3; i++)
    {
        (*k)[i] = (__float128*) malloc (s * sizeof (__float128));
        for (j = 0; j < s; j++) (*k)[i][j] = 0;
    }
}

// CAB matrix reading
void reading (FILE* inpf, __float128** cab, unsigned int s,
              long long unsigned i, long long unsigned j)
{
    for (i = 0; i < s; i++) cab[0][i] = rd (inpf);

    for (i = 2; i < s + 3; i++)
        for (j = 0; j + 1 < i && j < s && !feof (inpf); j++)
            cab[i][j] = rd (inpf);
}

// Reading floating point number from a/b format
__float128 rd (FILE* inpf)
{
    long long int chisl = 0;
    long long int znam  = 1;
    char          tmp;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"

    fscanf (inpf, "%Ld", &chisl);
    if (fscanf (inpf, "%c", &tmp) && tmp == '/') fscanf (inpf, "%Ld", &znam);

#pragma GCC diagnostic pop

    return ((__float128) chisl) / znam;
}

// CAB matrix printing
void printing (__float128** cab, unsigned int s, long long unsigned i,
               long long unsigned j)
{
    printf ("\n");
    for (i = 0; i < s + 3; i++)
    {
        for (j = 0; j < s; j++) printf ("%6.3Lf ", (long double) cab[i][j]);
        printf ("\n");
    }
    printf ("\n");
}

// Deconstructor
void cleaning (long long unsigned i, unsigned int s, FILE* inpf,
               __float128** cab, __float128** k)
{
    fclose (inpf);

    for (i = 0; i < s + 3; i++) free (cab[i]);
    free (cab);

    for (i = 0; i < 3; i++) free (k[i]);
    free (k);
}
