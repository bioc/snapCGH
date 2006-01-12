#ifdef MSDOS		/* 80x87 math chip microsoft C 	*/
#define EPSILON 1.0e-12
#define SQREPSILON 1.0e-6
#else			/* vax 11/750 unix cc		*/
#define EPSILON 1.0e-12
#define SQREPSILON 1.0e-6
#endif

/* dimension of arrays - can be increased*/
#define N 35

typedef struct
{
  double *data;
  double *covars1;
  int nrow;
}
Exts;
