
//Structure to pass data to the functions to be evaluated.

typedef struct
{
  double *data;
  double *covars1;
  double *covars2;
  double *covars3;
  int ncovars;
  int var;
  int nrow;
}
dataStore;  

//Function prototypes

void runNelderMead(int*, double*, double*, double*, double*, double*, double*, double*, int*, int*, double*, int*, int*, int*);
  static double fr_five(int, double[35], void*);
  static double fr_four(int, double[24], void*);
  static double fr_three(int, double[15], void*);
  static double fr_two(int, double[8], void*);

