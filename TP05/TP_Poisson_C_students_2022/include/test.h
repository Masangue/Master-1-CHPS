#ifndef TEST_H
#define TEST_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "atlas_headers.h"
#include "lib_poisson1D.h"

/*
    Test si l'implementation du format de stockage de la matrice est correcte

    Comapre le résultat entre dgemv et dgbmv
*/
int test_GB_storage_poisson1D(double *AB, int *la, int *lab, int *kl, int *ku, int *kv, double *X);

#endif