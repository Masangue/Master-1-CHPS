/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

/*
    Remplir le tableau 1D correspondant à la matrice poisson1D en general band, priorité colonne
*/
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv)
{
    int col;
    for (int ii=0;ii<(*la);ii++)
    {
        col = ii*(*lab);
        for (int jj=0;jj<(*kv);jj++)
        {
            AB[col + jj] = 0.0;
        }
        AB[col + (*kv)]     = -1.0;
        AB[col + (*kv) + 1] = 2.0;
        AB[col + (*kv) + 2] = -1.0;
    }
    AB[*kv] = 0.0;
    AB[(*lab)*(*la)-1] = 0.0;
}

/*
    Remplir le tableau 1D correspondant à la matrice identite en general band, priorité colonne
*/
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv)
{
    int col;
    for (int ii=0;ii<(*la);ii++)
    {
        col = ii*(*lab);
        for (int jj=0;jj<(*kv);jj++)
        {
            AB[col + jj] = 0.0;
        }
        AB[col + (*kv)]     = 0.0;
        AB[col + (*kv) + 1] = 1.0;
        AB[col + (*kv) + 2] = 0.0;
    }
}

/*
    Initialiser le vecteur RHS à 0 sauf aux extrémités qui valent les conditions aux limites BC0 et BC1
*/
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1)
{
    RHS[0] = *BC0;
    RHS[(*la)-1] = *BC1;

    for (int ii=1;ii<(*la)-1;ii++)
    {
        RHS[ii] = 0.0;
    }
}  

/*
    Le sujet donne comme solution analytique T(x) = T0 + x(T1-T0)
*/
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1)
{
    
    double f_diff = (*BC1)-(*BC0);
  
    for (int ii=0;ii<(*la);ii++)
    {
        EX_SOL[ii] = (*BC0) + f_diff*X[ii];
    }
}  

/*
    Initialiser le vecteur X avec les points de discrétisation
*/
void set_grid_points_1D(double* x, int* la)
{
  int n = (*la)+1;

  for (int ii=0;ii<(*la);ii++)
  {
      x[ii] = (ii+1.0)/n;
  }
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
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

int indexABCol(int i, int j, int *lab){
  return i*(*lab)+j;
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