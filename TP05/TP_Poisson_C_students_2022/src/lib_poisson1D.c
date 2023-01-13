/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"

/*
    Remplir le tableau 1D correspondant à la matrice poisson1D en general band, priorité colonne
*/
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}

/*
    Remplir le tableau 1D correspondant à la matrice identite en general band, priorité colonne
*/
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

/*
    Initialiser le vecteur RHS à 0 sauf aux extrémités qui valent les conditions aux limites BC0 et BC1
*/
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

/*
    Le sujet donne comme solution analytique T(x) = T0 + x(T1-T0)
*/
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

/*
    Initialiser le vecteur X avec les points de discrétisation
*/
void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void eig_poisson1D(double* eigval, int *la){
}

double eigmax_poisson1D(int *la){
  return 0;
}

double eigmin_poisson1D(int *la){
  return 0;
}

double richardson_alpha_opt(int *la){
  return 0;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){

}

/*
    Extract the D matrix from the AB matrix for Jacobi method
*/
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv)
{
    // extract the diagonal matrix D from AB
    for (int ii=0;ii<(*la);ii++)
    { 
        for (int jj=0;jj<(*kv);jj++)
        {
            MB[ii * (*lab) + jj] = 0.0;
        }
        MB[ii * (*lab) + (*kv)    ] = 0.0;
        MB[ii * (*lab) + (*kv) + 1] = AB[ii * (*lab) + (*kv) + 1];
        MB[ii * (*lab) + (*kv) + 2] = 0.0;
    }
}

/*
    Extract the D-E matrix from the AB matrix for Gauss-Seidel method
*/
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv)
{
    // extract the lower triangular matrix D-E from AB
    for (int ii=0;ii<(*la);ii++)
    { 
        for (int jj=0;jj<(*kv);jj++)
        {
            MB[ii * (*lab) + jj] = 0.0;
        }
        MB[ii * (*lab) + (*kv)    ] = 0.0;
        MB[ii * (*lab) + (*kv) + 1] = AB[ii * (*lab) + (*kv) + 1];  // D matrix
        MB[ii * (*lab) + (*kv) + 2] = -AB[ii * (*lab) + (*kv) + 2]; // -E matrix
    }
}


/*
    Algoritme itératif par Méthode de Richardson

    @param resvec : vecteur des résidus relatifs à chaque itération
*/
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite)
{
    double *temp_vec = (double*)malloc((*la)*sizeof(double));
    double norm_RHS = cblas_dnrm2((*la), RHS, 1); 
    double norm_bmAx = 0.0;

    // factorise the matrix MB in LU for inversion
    int info = 0;
    int nrhs = 1;
    int ku_local = 0;

    int *ipiv = (int*)malloc((*la)*sizeof(int));
    double *MB_LU = (double*)malloc((*la)*(*lab)*sizeof(double));
    cblas_dcopy((*la)*(*lab), MB, 1, MB_LU, 1);

    dgbtrf_(la, la, kl, &ku_local, MB_LU, lab, ipiv, &info);

    // compute the norm of the residual at step 0
    cblas_dgbmv(CblasColMajor, CblasNoTrans, (*la), (*la), (*kl), (*ku), 1.0, AB, (*lab), X, 1, 0.0, temp_vec, 1);  // compute A.X_0
    cblas_daxpy((*la), -1.0, RHS, 1, temp_vec, 1); // compute -b+A.X_0 result
    resvec[0] = cblas_dnrm2((*la), temp_vec, 1);   // compute norm of b+A.X_0 (=norm -b+A.X_0)

    for (int kk=1;(kk<(*maxit) && resvec[kk-1]>(*tol));kk++)
    {
        // compute the solution at step kk
        dgbtrs_("N", la, kl, &ku_local, &nrhs, MB_LU, lab, ipiv, temp_vec, la, &info, 1); // compute M^-1.(b+A.X_kk)
        cblas_daxpy((*la), -1.0, temp_vec, 1, X, 1); // compute X_kk+1 = X_kk - M^-1.(b+A.X_kk)
        
        // compute the norm of the residual at step kk
        cblas_dgbmv(CblasColMajor, CblasNoTrans, (*la), (*la), (*kl), (*ku), 1.0, AB, (*lab), X, 1, 0.0, temp_vec, 1);  // compute A.X_kk
        cblas_daxpy((*la), -1.0, RHS, 1, temp_vec, 1); // compute -b+A.X_kk result
        norm_bmAx = cblas_dnrm2((*la), temp_vec, 1);   // compute norm of b+A.X_kk (=norm -b+A.X_kk)
        resvec[kk] =  norm_bmAx / norm_RHS;            // compute relative residual
    }
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

/*
    Factorisation LU de la matrice tridiagonale AB
    @param ipiv : vecteur des pivots, mit à ipiv[i] = i+1 pour dgbtrs
*/
int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info)
{
    ipiv[0] = 1;
    for (int ii=1;ii<(*n);ii++)
    {
      AB[ii*(*lab)+2] -= (AB[(ii-1)*(*lab)+3]*AB[ii*(*lab)+1]) / AB[(ii-1)*(*lab)+2];
      AB[(ii-1)*(*lab)+3] /= AB[(ii-1)*(*lab)+2];
      ipiv[ii] = ii+1;
    }
}
