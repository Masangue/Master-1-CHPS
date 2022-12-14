#######################################
# masangue-victus-16-d0195nf.mk
# Default options for masangue computer
#######################################
CC=gcc
LIBSLOCAL=-llapack -lcblas -lm
INCLUDEBLASLOCAL=-I/usr/include
OPTCLOCAL=-fPIC -march=native