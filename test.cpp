#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//Function for MOM estimator
void test(int* nrow, int* ncol, double* dat, int* nrow_Dmat, double* sDmat, int* length_Hvec, double* Hvec, double* eps, int* t, double* G, double* N){

  double** mdat = (double**) malloc(sizeof(double*) + *ncol);
  for(int i = 0; i < *ncol; i++) mdat[i] = dat + i * *nrow;

  double** mDmat = (double**) malloc(sizeof(double*) + *nrow_Dmat);
  for(int i = 0; i < *nrow_Dmat; i++) mDmat[i] = sDmat + i * *nrow_Dmat;

  double** mG = (double**) malloc(sizeof(double*) + *length_Hvec);
  for(int i = 0; i < *length_Hvec; i++) mG[i] = G + i * *t;

  double** mN = (double**) malloc(sizeof(double*) + *length_Hvec);
  for(int i = 0; i < *length_Hvec; i++) mN[i] = N + i * *t;


    for(int i = 0; i < *nrow; i++){
      for(int j = i; j < *nrow; j++){

       //spherical distance of a pair
       double sDist = mDmat[(int) mdat[i][4]][(int) mdat[j][4]];

       //time lag of a pair
       int tDist = (int) abs(mdat[i][2] - mdat[j][2]);

       for(int k = 0; k < *t; k++){
        for(int l = 0; l < *length_Hvec; l++){
          //if((Hvec[l] - *eps) < sDist & sDist <= (Hvec[l] + *eps) & tDist == k){
            //#number of obs
            mN[l][k] += 1;
            //#sum of covariances to compute MOM estimate
            mG[l][k] += mdat[i][3] * mdat[j][3];
        //  }
        }
       }
      }
    }


  free(mdat);
  free(mDmat);
  free(mG);
  free(mN);
}



int main(){

  double data2[15][5] = {
     {2.465278, -1.76714587, 1, -3.78, 0},
     {1.854412, -0.06544985, 1,  0.64, 1},
     {2.596177, -1.24354709, 1,  4.98, 2},
     {1.199914,  1.67987941, 1,  0.27, 3},
     {2.465278, -2.72707696, 1,  0.43, 4},
     {2.465278, -1.76714587, 2,  0.21, 0},
     {1.854412, -0.06544985, 2, -0.82, 1},
     {2.596177, -1.24354709, 2,  0.44, 2},
     {1.199914,  1.67987941, 2, -0.01, 3},
     {2.465278, -2.72707696, 2, -1.11, 4},
     {2.465278, -1.76714587, 3, -3.84, 0},
     {1.854412, -0.06544985, 3, -0.02, 1},
     {2.596177, -1.24354709, 3,  1.28, 2},
     {1.199914,  1.67987941, 3,  0.65, 3},
     {2.465278, -2.72707696, 3,  0.02, 4}
 };

/*
 double dist[5][5] = {
   {0.2, 0.3, 0.5, 2.3, 3.7},
   {0.2, 0.3, 0.5, 2.3, 3.7},
   {0.2, 0.3, 0.5, 2.3, 3.7},
   {0.2, 0.3, 0.5, 2.3, 3.7},
   {0.2, 0.3, 0.5, 2.3, 3.7}
 };
 */

double dist2[5][5] = {
  {0.0000000, 9.613270e-01, 4.461393e-01, 5.491988e-01, 0.7367728},
  {0.9613270, 1.490116e-08, 9.402363e-01, 4.217217e-01, 1.1512551},
  {0.4461393, 9.402363e-01, 1.490116e-08, 6.664704e-01, 1.1779473},
  {0.5491988, 4.217217e-01, 6.664704e-01, 1.490116e-08, 0.8387904},
  {0.7367728, 1.151255e+00, 1.177947e+00, 8.387904e-01, 0.0000000}
};


 double Nmat2[10][3]={
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0}
 };

 double Gmat2[10][3]={
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0},
   {0, 0, 0}
 };



   double data[] = {2.46527757,  1.85441233,  2.59617726,  1.19991386,  2.46527757,  2.46527757,  1.85441233,  2.59617726,  1.19991386,
 2.46527757,  2.46527757,  1.85441233,  2.59617726,  1.19991386,  2.46527757, -1.76714587, -0.06544985, -1.24354709,
 1.67987941, -2.72707696, -1.76714587, -0.06544985, -1.24354709,  1.67987941, -2.72707696, -1.76714587, -0.06544985,
 -1.24354709,  1.67987941, -2.72707696,  1.00000000,  1.00000000,  1.00000000,  1.00000000,  1.00000000,  2.00000000,
 2.00000000,  2.00000000,  2.00000000,  2.00000000,  3.00000000,  3.00000000,  3.00000000,  3.00000000,  3.00000000,
 -3.78000000,  0.64000000,  4.98000000,  0.27000000,  0.43000000,  0.21000000, -0.82000000,  0.44000000, -0.01000000,
 -1.11000000, -3.84000000, -0.02000000,  1.28000000,  0.65000000,  0.02000000,  0.00000000,  1.00000000,  2.00000000,
 3.00000000,  4.00000000,  0.00000000,  1.00000000,  2.00000000,  3.00000000,  4.00000000,  0.00000000,  1.00000000,
 2.00000000,  3.00000000,  4.00000000};


 double dist[] = {0.000000e+00, 9.613270e-01, 4.461393e-01, 5.491988e-01, 7.367728e-01, 9.613270e-01, 1.490116e-08, 9.402363e-01,
   4.217217e-01, 1.151255e+00, 4.461393e-01, 9.402363e-01, 1.490116e-08, 6.664704e-01, 1.177947e+00, 5.491988e-01,
   4.217217e-01, 6.664704e-01, 1.490116e-08, 8.387904e-01, 7.367728e-01, 1.151255e+00, 1.177947e+00, 8.387904e-01,
   0.000000e+00};

 double Nmat[30]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
 double Gmat[30]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

 double Hvector[10] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
 double epsilon = 0.1;

 double sDist=0;
 int tDist=0;



 int nrow1 = 15, ncol1 = 5, nrow2=5, len_Hvec=10, tt=3;
 int *nrow = &nrow1, *ncol = &ncol1, *nrow_Dmat=&nrow2, *length_Hvec=&len_Hvec, *t=&tt;
 double *Hvec = &Hvec[0], *eps=&epsilon;
 double *dat = &data[0], *sDmat = &dist[0], *G = & Gmat[0], *N = &Nmat[0];

 double** mdat = (double**) malloc(sizeof(double*) * *ncol);
 for(int i = 0; i < *ncol; i++) mdat[i] = dat + i * *nrow;

 double** mDmat = (double**) malloc(sizeof(double*) * *nrow_Dmat);
 for(int i = 0; i < *nrow_Dmat; i++) mDmat[i] = sDmat + i * *nrow_Dmat;

 double** mG = (double**) malloc(sizeof(double*) * *length_Hvec);
 for(int i = 0; i < *length_Hvec; i++) mG[i] = G + i * *t;

 double** mN = (double**) malloc(sizeof(double*) * *length_Hvec);
 for(int i = 0; i < *length_Hvec; i++) mN[i] = N + i * *t;



 for(int i = 0; i < *nrow; i++){
   for(int j = i; j < *nrow; j++){

    //spherical distance of a pair
    sDist = mDmat[(int) mdat[4][i]][(int) mdat[4][j]];
    //std::cout << (int) mdat[i][4] << "\n";
    //std::cout << (int) mdat[j][4] << "\n";
    //std::cout << sDist << "\n\n";

    //time lag of a pair
    tDist = (int) abs(mdat[2][i] - mdat[2][j]);
    //std::cout << tDist << "\n";

    for(int k = 0; k < *t; k++){
     for(int l = 0; l < *length_Hvec; l++){
       if((Hvec[l] - *eps) < sDist & sDist <= (Hvec[l] + *eps) & tDist == k){
         //#number of obs
         mN[k][l] += 1;
         //#sum of covariances to compute MOM estimate
         mG[k][l] += mdat[3][i] * mdat[3][i];

         //std::cout << mG[l][k] << "\n";
         //std::cout << mN[l][k] << "\n";
       }
     }
    }
   }
 }

 /*
  std::cout << "1 and 2" << "\n";
  std::cout << data1[1][2] << "\n";
  std::cout << mdat[2][1] << "\n";
  */
  std::cout << "0 and 1" << "\n";
  std::cout << data2[13][2] << "\n";
  std::cout << mdat[2][13] << "\n\n";

  std::cout << mN[1][2] << "\n";
  std::cout << mG[1][2] << "\n";


 //std::cout << (int) mdat[0][4] << "\n";
 //std::cout << mDmat[(int) mdat[14][4]][(int) mdat[0][4]] << "\n";

 free(mdat);
 free(mDmat);
 free(mG);
 free(mN);

}
