##########################################
# Makefile                               #
# Makefile for the code developed in TP1 #
#                                        #
# T. Dufaud                              #
##########################################
################################
# Variables for this makefile
################################
# 
# -- option to dedicated machine
#
# Rk: You should create a file such as your-machineName.mk
# Follow the example of the machine called "ambre" in the 
# file ambre.mk
#
HOSTNAME=$(shell hostname)
include $(HOSTNAME).mk

# 
# -- Compiler Option
OPTC=${OPTCLOCAL}

#
# -- Directories
TPDIR=.
TPDIRSRC=$(TPDIR)/src

#
# -- librairies
LIBS=${LIBSLOCAL}

# -- Include directories
INCLATLAS=${INCLUDEBLASLOCAL}
INCL= -I $(TPDIR)/include $(INCLATLAS)

#
#################################################################
# makefile
############
#
OBJENV= obj/tp_env.o
OBJTP2ITER= obj/lib_poisson1D.o obj/tp_poisson1D_iter.o
OBJTP2DIRECT= obj/lib_poisson1D.o obj/tp_poisson1D_direct.o
#

all: bin/tp_testenv bin/tpPoisson1D_iter bin/tpPoisson1D_direct

testenv: bin/tp_testenv

tpPoisson1D_iter: bin/tpPoisson1D_iter

tpPoisson1D_direct: bin/tpPoisson1D_direct

obj/tp_env.o: $(TPDIRSRC)/tp_env.c
	@mkdir -p obj
	$(CC) $(OPTC) -o $@ -c $(INCL) $(TPDIRSRC)/tp_env.c 

obj/lib_poisson1D.o: $(TPDIRSRC)/lib_poisson1D.c
	@mkdir -p obj
	$(CC) $(OPTC) -o $@ -c $(INCL) $(TPDIRSRC)/lib_poisson1D.c 

obj/tp_poisson1D_iter.o: $(TPDIRSRC)/tp_poisson1D_iter.c
	@mkdir -p obj
	$(CC) $(OPTC) -o $@ -c $(INCL) $(TPDIRSRC)/tp_poisson1D_iter.c  

obj/tp_poisson1D_direct.o: $(TPDIRSRC)/tp_poisson1D_direct.c
	@mkdir -p obj
	$(CC) $(OPTC) -o $@ -c $(INCL) $(TPDIRSRC)/tp_poisson1D_direct.c  

bin/tp_testenv: $(OBJENV)
	@mkdir -p bin
	$(CC) -o bin/tp_testenv $(OPTC) $(OBJENV) $(LIBS)

bin/tpPoisson1D_iter: $(OBJTP2ITER)
	@mkdir -p bin
	$(CC) -o bin/tpPoisson1D_iter $(OPTC) $(OBJTP2ITER) $(LIBS)

bin/tpPoisson1D_direct: $(OBJTP2DIRECT)
	@mkdir -p bin
	$(CC) -o bin/tpPoisson1D_direct $(OPTC) $(OBJTP2DIRECT) $(LIBS)

run_testenv:
	bin/tp_testenv

run_tpPoisson1D_iter:
	bin/tpPoisson1D_iter

run_tpPoisson1D_direct:
	bin/tpPoisson1D_direct

clean:
	rm *.o bin/*
