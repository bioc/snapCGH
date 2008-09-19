#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <R_ext/Applic.h>
#include "header.h"

static double fr_two (int n, double par[8], void *ex){

  dataStore *ext;
  double *data, *covars1, *covars2, *covars3;
  int j, k, t, i;
  ext = ex;

  unsigned int nrow1 = ext->nrow;
  unsigned int ncovars = ext->ncovars;
  unsigned int varfixed = ext->var;
  data = ext->data;
  covars1 = ext->covars1;
  covars2 = ext->covars2;
  covars3 = ext->covars3; 

  double pr1[2], p1, p2, rate1;
  double S[2], mu[2], Sigma[2];
  double gammaA[2][2], gammaB[2][2], gammaC[2][2], alpha[2][nrow1], emis_prob[2][nrow1];
  double denom, temp, temp2, temp3;

  //initialize output
  double output = 0;

  mu[0] = par[0];
  mu[1] = par[1];
  Sigma[0] = par[2];
  Sigma[1] = par[3];
  double prior = par[4];
  double eta = par[5];
  double zeta = par[6];
  double omega = par[7];
  S[0] = exp(Sigma[0]);

  if(varfixed == 0){
    S[1] = exp(Sigma[1]);
  }
  else{
    S[1] = S[0];
      }
  
  if (prior < 20) {
    pr1[0] = exp(prior)/(1+exp(prior)); 
      }
  else {
    pr1[0] = 1;
      }

  pr1[1] = 1 - pr1[0];
  
  if (eta < 20) {
    p1 = exp(eta)/(1+exp(eta));
      }
  else {
    p1 = 1;
      }

  if (zeta < 20) {
    p2 = exp(zeta)/(1+exp(zeta)); 
      }
  else {
    p2 = 1;
      }

  if (omega < 20) {
    rate1 = exp(omega); 
      }
  else {
    rate1 = exp(20);
      }

  gammaA[0][0] = (1-p1); gammaA[0][1] = p1; gammaA[1][0] = p2; gammaA[1][1] = (1-p2);
  gammaB[0][0] = p1; gammaB[0][1] = -p1; gammaB[1][0] = -p2; gammaB[1][1] = p2;
  
  for(j = 0; j < nrow1; j++){
    //for(k = 0; k < 2; k++){
    for(k = 2; k--;){
      if(S[k] > 0.001){
	emis_prob[k][j] = dnorm(data[j], mu[k], S[k], 0);
      }
      else {
	if(data[j] >= mu[k]*(0.999) && data[j] <= mu[k]*(1.001)){
	  emis_prob[k][j] = 1;
	}
	else {
	  emis_prob[k][j] = 0;
	}
      }
    }
  }

  for(t = 1; t < nrow1; t++){

    switch(ncovars){
    case 1:
      temp = (exp(-(pow(covars1[t-1],rate1))));
      gammaC[0][0] = gammaA[0][0] + (temp*gammaB[0][0]);
      gammaC[0][1] = gammaA[0][1] + (temp*gammaB[0][1]);
      gammaC[1][0] = gammaA[1][0] + (temp*gammaB[1][0]);
      gammaC[1][1] = gammaA[1][1] + (temp*gammaB[1][1]);
      break;
    case 2:
      temp = exp(-(pow(covars1[t-1],rate1))*covars2[t-1]);
      gammaC[0][0] = gammaA[0][0] + (temp*gammaB[0][0]);
      gammaC[0][1] = gammaA[0][1] + (temp*gammaB[0][1]);
      gammaC[1][0] = gammaA[1][0] + (temp*gammaB[1][0]);
      gammaC[1][1] = gammaA[1][1] + (temp*gammaB[1][1]);
      break;
    case 3:
      temp = exp(-(pow(covars1[t-1],rate1))*covars2[t-1]*covars3[t-1]);
      gammaC[0][0] = gammaA[0][0] + (temp*gammaB[0][0]);
      gammaC[0][1] = gammaA[0][1] + (temp*gammaB[0][1]);
      gammaC[1][0] = gammaA[1][0] + (temp*gammaB[1][0]);
      gammaC[1][1] = gammaA[1][1] + (temp*gammaB[1][1]);
      break;
    }

    // compute alpha[aa][0]
    
    double c[nrow1];
    
    c[0] = 0;
    for (i=0;i<2;i++){
        alpha[i][0] = pr1[i]*emis_prob[i][0];
        c[0] = c[0] + alpha[i][0]; 
    }

    // scale the alpha[i][0]

    c[0] = 1/c[0];
    for (j=0;j<2;j++){
        alpha[j][0] = c[0] * alpha[j][0];
    }

    // compute alpha[i][t] (t = exon, i = state) 
    
    int d,l;
    
    for(t = 1; t < nrow1; t++){
        c[t] = 0;
        for(k = 0; k < 2; k++){
            alpha[k][t] = 0;
            for(d=0;d<2;d++){
                alpha[k][t] = alpha[k][t] + alpha[d][t-1]*gammaC[d][k];
            }
            alpha[k][t] = alpha[k][t]*emis_prob[k][t];
            c[t] = c[t] + alpha[k][t];
        }

        // scale alpha[k][t]

        c[t] = 1/c[t];
        
        for (l=0;l<2;l++){
            alpha[l][t] = c[t]*alpha[l][t];
        }
    }
        
    // log-likelihood

    int q;
    output = 0;
    for (q=0;q<nrow1;q++){
        output = output + log(c[q]);
    }
    output = ((-1)*output);
  }    
    return(-1*(output));
}
    
static double fr_three(int n, double par[15], void *ex){

  dataStore *ext;
  double *data, *covars1, *covars2, *covars3;
  int j, k, t, i, m;
  ext = ex;

  unsigned int nrow1 = ext->nrow;
  data = ext->data;
  unsigned int ncovars = ext->ncovars;
  covars1 = ext->covars1;
  covars2 = ext->covars2;
  covars3 = ext->covars3;
  unsigned int varfixed = ext->var;

  double pr[3], p1, p2, p3, p4, p5, p6, rate1;
  double S[3], mu[3], Sigma[3];
  double gammaA[3][3], gammaB[3][3], gammaC[3][3], alpha[3][nrow1], alphahat[3][nrow1], emis_prob[3][nrow1];
  double denom, temp, temp2, temp3;

  //initialize output
  double output = 0;
    
  mu[0] = par[0];
  mu[1] = par[1];
  mu[2] = par[2];
  Sigma[0] = par[3];
  Sigma[1] = par[4];
  Sigma[2] = par[5];
  double prior1 = par[6];
  double prior2 = par[7];
  double eta = par[8];
  double zeta = par[9];
  double theta = par[10];
  double beta = par[11];
  double gamma = par[12];
  double xi = par[13];
  double omega = par[14];
  S[0] = exp(Sigma[0]);
  if(varfixed == 0){    
    S[1] = exp(Sigma[1]);
    S[2] = exp(Sigma[2]);
  }
  else{
    S[1] = S[2] = S[0];
  }
  
  if (prior1 < 20) {
    pr[0] = exp(prior1)/(1+exp(prior1)); 
      }
  else {
    pr[0] = 1;
      }

   if (prior2 < 20) {
    pr[1] = (1 - pr[0])*exp(prior2)/(1+exp(prior2)); 
      }
  else {
    pr[1] = (1- pr[0]);
      }

   pr[2] = 1 - pr[0] - pr[1];
  
  if (eta < 20) {
    p1 = exp(eta)/(1+exp(eta));
      }
  else {
    p1 = 1;
      }

  if (zeta < 20) {
    p2 = (1 - p1)*exp(zeta)/(1+exp(zeta)); 
      }
  else {
    p2 = (1 - p1);
      }

  if (theta < 20) {
    p3 = exp(theta)/(1+exp(theta));
      }
  else {
    p3 = 1;
      }

  if (beta < 20) {
    p4 = (1 - p3)*exp(beta)/(1+exp(beta)); 
      }
  else {
    p4 = (1 - p3);
      }

  if (gamma < 20) {
    p5 = exp(gamma)/(1+exp(gamma));
      }
  else {
    p5 = 1;
      }

  if (xi < 20) {
    p6 = (1 - p5)*exp(xi)/(1+exp(xi)); 
      }
  else {
    p6 = (1 - p5);
      }

  if (omega < 20) {
    rate1 = exp(omega); 
      }
  else {
    rate1 = exp(20);
      }

  gammaA[0][0] = (1-p1-p2); gammaA[0][1] = p1; gammaA[0][2] = p2; gammaA[1][0] = p3; gammaA[1][1] = (1-p3-p4); gammaA[1][2] = p4; gammaA[2][0] = p5; gammaA[2][1] = p6; gammaA[2][2] = (1-p5-p6);
  gammaB[0][0] = (p1+p2); gammaB[0][1] = -p1; gammaB[0][2] = -p2; gammaB[1][0] = -p3; gammaB[1][1] = (p3+p4); gammaB[1][2] = -p4; gammaB[2][0] = -p5; gammaB[2][1] = -p6; gammaB[2][2] = (p5+p6);
  
  //  Rprintf("%f \t %f\n%f \t %f\n", gammaA[0][0],gammaA[0][1],gammaA[1][0],gammaA[1][1]);

  for(j = 0; j < nrow1; j++){
    for(k = 3; k--;){
      if(S[k] > 0.001){
	emis_prob[k][j] = dnorm(data[j], mu[k], S[k], 0);
      }
      else {
	if(data[j] >= mu[k]*(0.999) && data[j] <= mu[k]*(1.001)){
	  emis_prob[k][j] = 1;
	}
	else {
	  emis_prob[k][j] = 0;
	}
      }
    }
  }

  for(t = 1; t < nrow1; t++){

    switch(ncovars){
    case 1:
      temp = exp(-(pow(covars1[t-1],rate1)));
      for(j = 0; j < 3; j++){
	for(m = 0; m < 3; m++){
	  gammaC[j][m] = gammaA[j][m] + (temp*gammaB[j][m]);
	}
      }
      break;
    case 2:
      temp = exp(-(pow(covars1[t-1],rate1))*covars2[t-1]);
      for(j = 0; j < 3; j++){
	for(m = 0; m < 3; m++){
	  gammaC[j][m] = gammaA[j][m] + (temp*gammaB[j][m]);
	}
      }
      break;
    case 3:
      temp = exp(-(pow(covars1[t-1],rate1))*covars2[t-1]*covars3[t-1]);
      for(j = 0; j < 3; j++){
	for(m = 0; m < 3; m++){
	  gammaC[j][m] = gammaA[j][m] + (temp*gammaB[j][m]);
	}
      }
      break;
    }


        // compute alpha[aa][0]
    
    double c[nrow1];
    
    c[0] = 0;
    for (i=0;i<3;i++){
        alpha[i][0] = pr[i]*emis_prob[i][0];
        c[0] = c[0] + alpha[i][0]; 
    }

    // scale the alpha[i][0]

    c[0] = 1/c[0];
    for (j=0;j<3;j++){
        alpha[j][0] = c[0] * alpha[j][0];
    }

    // compute alpha[i][t] (t = exon, i = state) 
    
    int d,l;
    
    for(t = 1; t < nrow1; t++){
        c[t] = 0;
        for(k = 0; k < 3; k++){
            alpha[k][t] = 0;
            for(d=0;d<3;d++){
                alpha[k][t] = alpha[k][t] + alpha[d][t-1]*gammaC[d][k];
            }
            alpha[k][t] = alpha[k][t]*emis_prob[k][t];
            c[t] = c[t] + alpha[k][t];
        }

        // scale alpha[k][t]

        c[t] = 1/c[t];
        
        for (l=0;l<3;l++){
            alpha[l][t] = c[t]*alpha[l][t];
        }
    }
        
    // log-likelihood

    int q;
    output = 0;
    for (q=0;q<nrow1;q++){
        output = output + log(c[q]);
    }
    output = ((-1)*output);
  }    
    return(-1*(output));
}    


static double fr_four(int n, double par[24], void *ex){

  dataStore *ext;
  ext = ex;

  double *data, *covars1, *covars2, *covars3;
  int j, k, t, i, m;

  int nrow1 = ext->nrow;
  data = ext->data;
  int ncovars = ext->ncovars;
  covars1 = ext->covars1;
  covars2 = ext->covars2;
  covars3 = ext->covars3;
  int varfixed = ext->var;

  double pr[4], p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, rate1;
  double S[4], mu[4], Sigma[4];
  double gammaA[4][4], gammaB[4][4], gammaC[4][4], alpha[4][nrow1], alphahat[4][nrow1], emis_prob[4][nrow1];
  double denom, temp, temp2, temp3;

  //initialize output
  double output = 0;
    
  mu[0] = par[0];
  mu[1] = par[1];
  mu[2] = par[2];
  mu[3] = par[3];
  Sigma[0] = par[4];
  Sigma[1] = par[5];
  Sigma[2] = par[6];
  Sigma[3] = par[7];
  double prior1 = par[8];
  double prior2 = par[9];
  double prior3 = par[10];
  double eta = par[11];
  double zeta = par[12];
  double nu = par[13];
  double theta = par[14];
  double beta = par[15];
  double phi = par[16];
  double gamma = par[17];
  double delt = par[18];
  double epsilon = par[19];
  double lambda = par[20];
  double rho = par[21];
  double xi = par[22];
  double omega = par[23];
  S[0] = exp(Sigma[0]);
  if(varfixed == 0){
    S[1] = exp(Sigma[1]);
    S[2] = exp(Sigma[2]);
    S[3] = exp(Sigma[3]);
  }
  else{
    S[1] = S[2] = S[3] = S[0];
  }
  
  if (prior1 < 20) {
    pr[0] = exp(prior1)/(1+exp(prior1)); 
      }
  else {
    pr[0] = 1;
      }

   if (prior2 < 20) {
    pr[1] = (1 - pr[0])*exp(prior2)/(1+exp(prior2)); 
      }
  else {
    pr[1] = (1- pr[0]);
      }

   if (prior3 < 20) {
    pr[2] = (1 -pr[0]-pr[1])*exp(prior3)/(1+exp(prior3)); 
      }
  else {
    pr[2] = (1-pr[0]-pr[1]);
      }

   pr[3] = 1 - pr[0] - pr[1] - pr[2];
  
  if (eta < 20) {
    p1 = exp(eta)/(1+exp(eta));
      }
  else {
    p1 = 1;
      }

  if (zeta < 20) {
    p2 = (1 - p1)*exp(zeta)/(1+exp(zeta)); 
      }
  else {
    p2 = (1 - p1);
      }

  if (nu < 20) {
    p3 = (1 -p1-p2)*exp(nu)/(1+exp(nu)); 
      }
  else {
    p3 = (1 - p1 - p2);
      }

  if (theta < 20) {
    p4 = exp(theta)/(1+exp(theta));
      }
  else {
    p4 = 1;
      }

  if (beta < 20) {
    p5 = (1 - p4)*exp(beta)/(1+exp(beta)); 
      }
  else {
    p5 = (1 - p4);
      }

  if (phi < 20) {
    p6 = (1 - p4 - p5)*exp(phi)/(1+exp(phi)); 
      }
  else {
    p6 = (1 - p4 - p5);
      }

  if (gamma < 20) {
    p7 = exp(gamma)/(1+exp(gamma));
      }
  else {
    p7 = 1;
      }

  if (delt < 20) {
    p8 = (1 - p7)*exp(delt)/(1+exp(delt)); 
      }
  else {
    p8 = (1 - p7);
      }

  if (epsilon < 20) {
    p9 = (1 - p7 - p8)*exp(epsilon)/(1+exp(epsilon)); 
      }
  else {
    p9 = (1 - p7 - p8);
      }

  if (lambda < 20) {
    p10 = exp(lambda)/(1+exp(lambda));
      }
  else {
    p10 = 1;
      }

  if (rho < 20) {
    p11 = (1 - p10)*exp(rho)/(1+exp(rho)); 
      }
  else {
    p11 = (1 - p10);
      }

  if (xi < 20) {
    p12 = (1 - p10 - p11)*exp(xi)/(1+exp(xi)); 
      }
  else {
    p12 = (1 - p10 - p11);
      }

  if (omega < 20) {
    rate1 = exp(omega); 
      }
  else {
    rate1 = exp(20);
      }

  gammaA[0][0] = (1-p1-p2-p3); gammaA[0][1] = p1; gammaA[0][2] = p2; gammaA[0][3] = p3; gammaA[1][0] = p4; gammaA[1][1] = (1-p4-p5-p6); gammaA[1][2] = p5; gammaA[1][3] = p6; gammaA[2][0] = p7; gammaA[2][1] = p8; gammaA[2][2] = (1-p7-p8-p9); gammaA[2][3] = p9; gammaA[3][0] = p10; gammaA[3][1] = p11; gammaA[3][2] = p12; gammaA[3][3] = (1-p10-p11-p12);
  gammaB[0][0] = (p1+p2+p3); gammaB[0][1] = -p1; gammaB[0][2] = -p2; gammaB[0][3] = -p3; gammaB[1][0] = -p4; gammaB[1][1] = (p4+p5+p6); gammaB[1][2] = -p5; gammaB[1][3] = -p6; gammaB[2][0] = -p7; gammaB[2][1] = -p8; gammaB[2][2] = (p7+p8+p9); gammaB[2][3] = -p9; gammaB[3][0] = -p10; gammaB[3][1] = -p11; gammaB[3][2] = -p12; gammaB[3][3] = (p10+p11+p12);

    for(j = 0; j < nrow1; j++){
     for(k = 4; k--;){
      if(S[k] > 0.001){
	emis_prob[k][j] = dnorm(data[j], mu[k], S[k], 0);
      }
      else {
	if(data[j] >= mu[k]*(0.999) && data[j] <= mu[k]*(1.001)){
	  emis_prob[k][j] = 1;
	}
	else {
	  emis_prob[k][j] = 0;
	}
      }
    }
  }

for(t = 1; t < nrow1; t++){  

    switch(ncovars){
    case 1:
      temp = exp(-(pow(covars1[t-1],rate1)));
      for(j = 0; j < 4; j++){
	for(m = 0; m < 4; m++){
	  gammaC[j][m] = gammaA[j][m] + (temp*gammaB[j][m]);
	}
      }
      break;
    case 2:
       temp = exp(-(pow(covars1[t-1],rate1))*covars2[t-1]);
      for(j = 0; j < 4; j++){
	for(m = 0; m < 4; m++){
	  gammaC[j][m] = gammaA[j][m] + (temp*gammaB[j][m]);
	}
      }
      break;
    case 3:
       temp = exp(-(pow(covars1[t-1],rate1))*covars2[t-1]*covars3[t-1]);
      for(j = 0; j < 4; j++){
	for(m = 0; m < 4; m++){
	  gammaC[j][m] = gammaA[j][m] + (temp*gammaB[j][m]);
	}
      }
      break;
    }


              // compute alpha[aa][0]
    
    double c[nrow1];
    
    c[0] = 0;
    for (i=0;i<4;i++){
        alpha[i][0] = pr[i]*emis_prob[i][0];
        c[0] = c[0] + alpha[i][0]; 
    }

    // scale the alpha[i][0]

    c[0] = 1/c[0];
    for (j=0;j<4;j++){
        alpha[j][0] = c[0] * alpha[j][0];
    }

    // compute alpha[i][t] (t = exon, i = state) 
    
    int d,l;
    
    for(t = 1; t < nrow1; t++){
        c[t] = 0;
        for(k = 0; k < 4; k++){
            alpha[k][t] = 0;
            for(d=0;d<4;d++){
                alpha[k][t] = alpha[k][t] + alpha[d][t-1]*gammaC[d][k];
            }
            alpha[k][t] = alpha[k][t]*emis_prob[k][t];
            c[t] = c[t] + alpha[k][t];
        }

        // scale alpha[k][t]

        c[t] = 1/c[t];
        
        for (l=0;l<4;l++){
            alpha[l][t] = c[t]*alpha[l][t];
        }
    }
        
    // log-likelihood

    int q;
    output = 0;
    for (q=0;q<nrow1;q++){
        output = output + log(c[q]);
    }
    output = ((-1)*output);
  }    
    return(-1*(output));
}


static double fr_five(int n, double par[35], void *ex){

  dataStore *ext;
  double *data, *covars1, *covars2, *covars3;
  int j, k, t, i, m;
  ext = ex;

  int nrow1 = ext->nrow;
  data = ext->data;
  int ncovars = ext->ncovars;
  covars1 = ext->covars1;
  covars2 = ext->covars2;
  covars3 = ext->covars3;
  int varfixed = ext->var;

  double  pr[5], p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, rate1;
  double S[5], mu[5], Sigma[5];
  double gammaA[5][5], gammaB[5][5], gammaC[5][5], alpha[5][nrow1], alphahat[5][nrow1], emis_prob[5][nrow1];
  double denom, temp, temp2, temp3;

  //initialize output
  double output = 0;
    
  mu[0] = par[0];
  mu[1] = par[1];
  mu[2] = par[2];
  mu[3] = par[3];
  mu[4] = par[4];
  Sigma[0] = par[5];
  Sigma[1] = par[6];
  Sigma[2] = par[7];
  Sigma[3] = par[8];
  Sigma[4] = par[9];
  double prior1 = par[10];
  double prior2 = par[11];
  double prior3 = par[12];
  double prior4 = par[13];
  double eta = par[14];
  double zeta = par[15];
  double nu = par[16];
  double omikron = par[17];
  double theta = par[18];
  double beta = par[19];
  double phi = par[20];
  double kappa = par[21];
  double gamma = par[22];
  double delt = par[23];
  double epsilon = par[24];
  double tau = par[25];
  double lambda = par[26];
  double rho = par[27];
  double xi = par[28];
  double iota = par[29];
  double chi = par[30];
  double upsilon = par[31];
  double psi = par[32];
  double aux = par[33];
  double omega = par[34];
  S[0] = exp(Sigma[0]);
  if(varfixed == 0){
    S[1] = exp(Sigma[1]);
    S[2] = exp(Sigma[2]);
    S[3] = exp(Sigma[3]);
    S[4] = exp(Sigma[4]);
  }
  else{
    S[1] = S[2] = S[3] = S[4] = S[0];
  }
  
  if (prior1 < 20) {
    pr[0] = exp(prior1)/(1+exp(prior1)); 
      }
  else {
    pr[0] = 1;
      }

   if (prior2 < 20) {
    pr[1] = (1 - pr[0])*exp(prior2)/(1+exp(prior2)); 
      }
  else {
    pr[1] = (1- pr[0]);
      }

   if (prior3 < 20) {
    pr[2] = (1 -pr[0]-pr[1])*exp(prior3)/(1+exp(prior3)); 
      }
  else {
    pr[2] = (1-pr[0]-pr[1]);
      }

   if (prior4 < 20) {
    pr[3] = (1 -pr[0]-pr[1]-pr[2])*exp(prior4)/(1+exp(prior4)); 
      }
  else {
    pr[3] = (1-pr[0]-pr[1]-pr[2]);
      }

   pr[4] = 1-pr[0]-pr[1]-pr[2]-pr[3];
   
  if (eta < 20) {
    p1 = exp(eta)/(1+exp(eta));
      }
  else {
    p1 = 1;
      }

  if (zeta < 20) {
    p2 = (1 - p1)*exp(zeta)/(1+exp(zeta)); 
      }
  else {
    p2 = (1 - p1);
      }

  if (nu < 20) {
    p3 = (1 -p1-p2)*exp(nu)/(1+exp(nu)); 
      }
  else {
    p3 = (1 - p1 - p2);
      }

  if (omikron < 20) {
    p4 = (1 -p1-p2-p3)*exp(omikron)/(1+exp(omikron)); 
      }
  else {
    p4 = (1 - p1 - p2 - p3);
      }

  if (theta < 20) {
    p5 = exp(theta)/(1+exp(theta));
      }
  else {
    p5 = 1;
      }

  if (beta < 20) {
    p6 = (1 - p5)*exp(beta)/(1+exp(beta)); 
      }
  else {
    p6 = (1 - p5);
      }

  if (phi < 20) {
    p7 = (1 - p5 - p6)*exp(phi)/(1+exp(phi)); 
      }
  else {
    p7 = (1 - p5 - p6);
      }

  if (kappa < 20) {
    p8 = (1 - p5 - p6 - p7)*exp(kappa)/(1+exp(kappa)); 
      }
  else {
    p8 = (1 - p5 - p6 - p7);
      }

  if (gamma < 20) {
    p9 = exp(gamma)/(1+exp(gamma));
      }
  else {
    p9 = 1;
      }

  if (delt < 20) {
    p10 = (1 - p9)*exp(delt)/(1+exp(delt)); 
      }
  else {
    p10 = (1 - p9);
      }

  if (epsilon < 20) {
    p11 = (1 - p9 - p10)*exp(epsilon)/(1+exp(epsilon)); 
      }
  else {
    p11 = (1 - p9 - p10);
      }

  if (tau < 20) {
    p12 = (1 - p9 - p10 - p11)*exp(tau)/(1+exp(tau)); 
      }
  else {
    p12 = (1 - p9 - p10 - p11);
      }

  if (lambda < 20) {
    p13 = exp(lambda)/(1+exp(lambda));
      }
  else {
    p13 = 1;
      }

  if (rho < 20) {
    p14 = (1 - p13)*exp(rho)/(1+exp(rho)); 
      }
  else {
    p14 = (1 - p13);
      }

  if (xi < 20) {
    p15 = (1 - p13 - p14)*exp(xi)/(1+exp(xi)); 
      }
  else {
    p15 = (1 - p13 - p14);
      }

  if (iota < 20) {
    p16 = (1 - p13 - p14 - p15)*exp(iota)/(1+exp(iota)); 
      }
  else {
    p16 = (1 - p13 - p14 - p15);
      }

  if (chi < 20) {
    p17 = exp(chi)/(1+exp(chi));
      }
  else {
    p17 = 1;
      }

  if (upsilon < 20) {
    p18 = (1 - p17)*exp(upsilon)/(1+exp(upsilon)); 
      }
  else {
    p18 = (1 - p17);
      }

  if (psi < 20) {
    p19 = (1 - p17 - p18)*exp(psi)/(1+exp(psi)); 
      }
  else {
    p19 = (1 - p17 - p18);
      }

  if (aux < 20) {
    p20 = (1 - p17 - p18 - p19)*exp(aux)/(1+exp(aux)); 
      }
  else {
    p20 = (1 - p17 - p18 - p19);
      }

  if (omega < 20) {
    rate1 = exp(omega); 
      }
  else {
    rate1 = exp(20);
      }

  gammaA[0][0] = (1-p1-p2-p3-p4); gammaA[0][1] = p1; gammaA[0][2] = p2; gammaA[0][3] = p3; gammaA[0][4] = p4;  
  gammaA[1][0] = p5; gammaA[1][1] = (1-p5-p6-p7-p8); gammaA[1][2] = p6; gammaA[1][3] = p7; gammaA[1][4] = p8; 
  gammaA[2][0] = p9; gammaA[2][1] = p10; gammaA[2][2] = (1-p9-p10-p11-p12); gammaA[2][3] = p11; gammaA[2][4] = p12; 
  gammaA[3][0] = p13; gammaA[3][1] = p14; gammaA[3][2] = p15; gammaA[3][3] = (1-p13-p14-p15-p16); gammaA[3][4] = p16; 
  gammaA[4][0] = p17; gammaA[4][1] = p18; gammaA[4][2] = p19; gammaA[4][3] = p20; gammaA[4][4] = (1-p17-p18-p19-p20);

  gammaB[0][0] = (p1+p2+p3+p4); gammaB[0][1] = -p1; gammaB[0][2] = -p2; gammaB[0][3] = -p3; gammaB[0][4] = -p4; 
  gammaB[1][0] = -p5; gammaB[1][1] = (p5+p6+p7+p8); gammaB[1][2] = -p6; gammaB[1][3] = -p7; gammaB[1][4] = -p8; 
  gammaB[2][0] = -p9; gammaB[2][1] = -p10; gammaB[2][2] = (p9+p10+p11+p12); gammaB[2][3] = -p11; gammaB[2][4] = -p12; 
  gammaB[3][0] = -p13; gammaB[3][1] = -p14; gammaB[3][2] = -p15; gammaB[3][3] = (p13+p14+p15+p16); gammaB[3][4] = -p16; 
  gammaB[4][0] = -p17; gammaB[4][1] = -p18; gammaB[4][2] = -p19; gammaB[4][3] = -p20; gammaB[4][4] = (p17+p18+p19+p20);

   for(j = 0; j < nrow1; j++){
    for(k = 5; k--;){
      if(S[k] > 0.001){
	emis_prob[k][j] = dnorm(data[j], mu[k], S[k], 0);
      }
      else {
	if(data[j] >= mu[k]*(0.999) && data[j] <= mu[k]*(1.001)){
	  emis_prob[k][j] = 1;
	}
	else {
	  emis_prob[k][j] = 0;
	}
      }
    }
  }

   for(t = 1; t < nrow1; t++){
  
    switch(ncovars){
    case 1:
      temp = exp(-(pow(covars1[t-1],rate1)));
      for(j = 0; j < 5; j++){
	for(m = 0; m < 5; m++){
	  gammaC[j][m] = gammaA[j][m] + (temp*gammaB[j][m]);
	}
      }
      break;
    case 2:
      temp = exp(-(pow(covars1[t-1],rate1))*covars2[t-1]);
      for(j = 0; j < 5; j++){
	for(m = 0; m < 5; m++){
	  gammaC[j][m] = gammaA[j][m] + (temp*gammaB[j][m]);
	}
      }
      break;
    case 3:
      temp = exp(-(pow(covars1[t-1],rate1))*covars2[t-1]*covars3[t-1]);
      for(j = 0; j < 5; j++){
	for(m = 0; m < 5; m++){
	  gammaC[j][m] = gammaA[j][m] + (temp*gammaB[j][m]);
	}
      }
      break;
    }
                  // compute alpha[aa][0]
    
    double c[nrow1];
    
    c[0] = 0;
    for (i=0;i<5;i++){
        alpha[i][0] = pr[i]*emis_prob[i][0];
        c[0] = c[0] + alpha[i][0]; 
    }

    // scale the alpha[i][0]

    c[0] = 1/c[0];
    for (j=0;j<5;j++){
        alpha[j][0] = c[0] * alpha[j][0];
    }

    // compute alpha[i][t] (t = exon, i = state) 
    
    int d,l;
    
    for(t = 1; t < nrow1; t++){
        c[t] = 0;
        for(k = 0; k < 5; k++){
            alpha[k][t] = 0;
            for(d=0;d<5;d++){
                alpha[k][t] = alpha[k][t] + alpha[d][t-1]*gammaC[d][k];
            }
            alpha[k][t] = alpha[k][t]*emis_prob[k][t];
            c[t] = c[t] + alpha[k][t];
        }

        // scale alpha[k][t]

        c[t] = 1/c[t];
        
        for (l=0;l<5;l++){
            alpha[l][t] = c[t]*alpha[l][t];
        }
    }
        
    // log-likelihood

    int q;
    output = 0;
    for (q=0;q<nrow1;q++){
        output = output + log(c[q]);
    }
    output = ((-1)*output);
  }    
    return(-1*(output));
 
}

void runNelderMead(int *nrow, double *xin, double *xout, double *Fmin, double *data, double *covars1, double *covars2, double *covars3, int *ncovars, int *var, double *epsilon, int *trace, int *numit, int *nstates){

  dataStore *ext;
  int fail;
  int fncount;

  ext = Calloc(1, dataStore);

  ext->data = data;
  ext->covars1 = covars1;
  ext->covars2 = covars2;
  ext->covars3 = covars3;
  ext->nrow = *nrow;
  ext->ncovars = *ncovars;
  
  double abstol = -HUGE_VAL;
  double intol = *epsilon;
  double alpha = 1;
  double beta = 0.5;
  double gamma = 2;

  switch(*nstates){
  case 2:
    nmmin(8, xin, xout, Fmin, fr_two, &fail, abstol, intol, ext, alpha, beta, gamma, *trace, &fncount, *numit);
    break;
  case 3:
    nmmin(15, xin, xout, Fmin, fr_three, &fail, abstol, intol, ext, alpha, beta, gamma, *trace, &fncount, *numit);
    break;
  case 4:
    nmmin(24, xin, xout, Fmin, fr_four, &fail, abstol, intol, ext, alpha, beta, gamma, *trace, &fncount, *numit);
    break;
  case 5:
    nmmin(35, xin, xout, Fmin, fr_five, &fail, abstol, intol, ext, alpha, beta, gamma, *trace, &fncount, *numit);
    break;
  default:
    Rprintf("Number of states not between 2 and 5\n");
    break;
  }

  free(ext);
}
