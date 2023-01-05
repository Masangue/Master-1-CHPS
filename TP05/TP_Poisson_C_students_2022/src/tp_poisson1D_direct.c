/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"
#include <string.h>

int main(int argc,char *argv[])
{
    int version = 1;

    if (argc == 2) {
        if (strcmp(argv[1], "dgbtrf") == 0) {
            version = 1;
        } else if (strcmp(argv[1], "dgbtrftridiag") == 0) {
            version = 2;
        } else if (strcmp(argv[1], "dgbsv") == 0) {
            version = 3;
        } else {
            printf("Usage: %s [dgbtrf|dgbtrftridiag|dgbsv]\n", argv[0]);
            exit(1);
        }
    }

    int ierr;
    int nbpoints, la;
    int ku, kl, kv, lab;
    int *ipiv;
    int info;
    int NRHS;
    double T0, T1;
    double *RHS, *EX_SOL, *X;
    double *AB;

    NRHS=1;
    nbpoints=10;
    la=nbpoints-2;
    T0=-5.0;
    T1=5.0;

    printf("--------- Poisson 1D ---------\n\n");
    RHS=(double *) malloc(sizeof(double)*la);
    EX_SOL=(double *) malloc(sizeof(double)*la);
    X=(double *) malloc(sizeof(double)*la);

    // TODO : you have to implement those functions
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
    
    write_vec(RHS, &la, "RHS.dat");
    write_vec(EX_SOL, &la, "EX_SOL.dat");
    write_vec(X, &la, "X_grid.dat");

    kv=1;
    ku=1;
    kl=1;
    lab=kv+kl+ku+1;

    AB = (double *) malloc(sizeof(double)*lab*la);

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

    // write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

    ipiv = (int *) calloc(la, sizeof(int));

    if (version == 1) 
    {
        printf("Solution with LAPACK dgbtrf + dgbtrs\n");

        /* LU Factorization */
        info=0;
        dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    //    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");

        /* Solution (Triangular) */
        if (info==0){
            dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info,1);
            if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
        }else{
            printf("\n INFO = %d\n",info);
        }
    } else if (version == 2) 
    {
        /* LU for tridiagonal matrix  (can replace dgbtrf_) */
        info=0;
        printf("Solution handmade for LU factorization of tridiagonal matrix\n");
        ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    //    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");

        /* Solution (Triangular) */
        if (info==0){
            dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info,1);
            if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
        }else{
            printf("\n INFO = %d\n",info);
        }

    } else if (version == 3) 
    {
        printf("Solution with LAPACK dgbsv\n");
        /* It can also be solved with dgbsv */
        dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    }
    
    // write_xy(RHS, X, &la, "SOL.dat");

    /* Relative forward error */
    // TODO : Compute relative norm of the residual
    double relres = 0.0;
    for (int i = 0; i < la; i++) {
        relres += fabs(RHS[i] - EX_SOL[i]) / fabs(EX_SOL[i]);
    }
    relres = relres / la;
    
    printf("\nThe relative forward error is relres = %e\n",relres);

    free(RHS);
    free(EX_SOL);
    free(X);
    free(AB);
    printf("\n\n--------- End -----------\n");
}
