#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "mom.h"

//Function for MOM estimator
void G_hat3(int* nrow, int* ncol, double* dat, int* nrow_Dmat, double* sDmat, int* length_Hvec, double* Hvec, double* eps, int* t, double* G, double* N){


  double sDist = 0;
  int tDist = 0;

  double** mdat = (double**) malloc(sizeof(double*) * *ncol);
  for(int i = 0; i < *ncol; i++) mdat[i] = dat + i * *nrow;

  double** mDmat = (double**) malloc(sizeof(double*) * *nrow_Dmat);
  for(int i = 0; i < *nrow_Dmat; i++) mDmat[i] = sDmat + i * *nrow_Dmat;

  double** mG = (double**) malloc(sizeof(double*) * *t);
  for(int i = 0; i < *t; i++) mG[i] = G + i * *length_Hvec;

  double** mN = (double**) malloc(sizeof(double*) * *t);
  for(int i = 0; i < *t; i++) mN[i] = N + i * *length_Hvec;


  for(int i = 0; i < *nrow; i++){
    for(int j = i; j < *nrow; j++){

     //spherical distance of a pair
     sDist = mDmat[(int) mdat[4][i]][(int) mdat[4][j]];

     //time lag of a pair
     tDist = (int) abs(mdat[2][i] - mdat[2][j]);

     for(int k = 0; k < *t; k++){
      for(int l = 0; l < *length_Hvec; l++){
        if((Hvec[l] - *eps) < sDist & sDist <= (Hvec[l] + *eps) & tDist == k){
          //#number of obs
          mN[k][l] += 1;
          //#sum of covariances to compute MOM estimate
          mG[k][l] += mdat[3][i] * mdat[3][j];
        }
      }
     }
    }
  }
  free(mdat);
  free(mDmat);
  free(mG);
  free(mN);
}
