#define EPSILON 1.0e-12
#define SQREPSILON 1.0e-6

/* dimension of arrays - can be increased*/
#define N 35

typedef struct
{
  double *data;
  double *covars1;
  int nrow;
}
Exts;
