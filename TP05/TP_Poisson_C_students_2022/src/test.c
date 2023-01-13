#include "test.h"

int test_GB_storage_poisson1D(double *AB, int *la, int *lab, int *kl, int *ku, int *kv, double *X)
{   
    // Compute result for AB in GB storage
    double * AB_res = (double *) malloc(sizeof(double)*(*la));

    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, (*ku)+(*kv), 1.0, AB, *lab, X, 1, 0.0, AB_res, 1);

    // A in dense format for Poisson 1D
    double *A = (double *) malloc(sizeof(double)*(*la)*(*la)); 
    for (int ii=0; ii<(*la)*(*la); ii++)
    {
       A[ii] = 0;
    }

    A[0] = 2;
    A[1] = -1;
    for (int ii=1; ii<(*la)-1; ii++)
    {   
        A[ii*(*la)+ii-1] = -1;
        A[ii*(*la)+ii] = 2;
        A[ii*(*la)+ii+1] = -1;
    }
    A[((*la)*(*la))-1] = 2;
    A[((*la)*(*la))-2] = -1;

    // Compute result for A in dense format
    double * A_res = (double *) malloc(sizeof(double)*(*la));

    cblas_dgemv(CblasColMajor, CblasNoTrans, *la, *la, 1.0, A, *la, X, 1, 0.0, A_res, 1);

    // Compare results
    int error = 0;
    for (int ii=0; ii<(*la); ii++)
    {
        if ( abs(AB_res[ii] - A_res[ii]) > 1e-15)
        {
            printf("## Error: AB_res[%d][%d] = %e != A_res[%d][%d] = %e\n", ii/(*la),ii%(*la), AB_res[ii], ii/(*la),ii%(*la), A_res[ii]);
            printf("## Something went wrong during test_GB_storage_poisson1D (dif : %lf)\n", AB_res[ii] - A_res[ii]);
            error = 1;
            break;
        }
    }

    free(AB_res);
    free(A_res);
    free(A);

    return error;
}


int main()
{

    int ierr;
    int jj;
    int nbpoints, la;
    int ku, kl, kv, lab;
    int *ipiv;
    int info;
    int NRHS;
    double T0, T1;
    double *RHS, *EX_SOL, *X;
    double **AAB;
    double *AB;

    double temp, relres;

    NRHS=1;
    nbpoints=10;
    la=nbpoints-2;
    T0=-5.0;
    T1=5.0;

    printf("--------------------Test---------------------\n");

    RHS=(double *) malloc(sizeof(double)*la);
    EX_SOL=(double *) malloc(sizeof(double)*la);
    X=(double *) malloc(sizeof(double)*la);

    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    kv=1;
    ku=1;
    kl=1;
    lab=kv+kl+ku+1;

    AB = (double *) malloc(sizeof(double)*lab*la);

    printf("## Testing GB storage for Poisson 1D ..... ##\n");

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

    double Xtest[la];
    srand(time(NULL));
    for (int i = 0; i < la; i++) Xtest[i] = rand() % 100;

    if (test_GB_storage_poisson1D(AB, &la, &lab, &kl, &ku, &ku, Xtest))
    {
        printf("## ................................ FAILED ##\n");
    } else {
        printf("## ................................ PASSED ##\n");
    }
}