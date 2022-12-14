all: clean test-env progtest run

# Check if cblas and lapack are installed
test-env:
	@echo CBLAS && if pacman -Q | grep -q cblas ; then echo found; else echo not found; fi;
	@echo ;
	@echo LAPACK && if pacman -Q | grep -q lapack ; then echo found; else echo not found; fi;
	@echo ;

# Try to compile with cblas and lapack to test if they are accessible
progtest: testenv.c
	gcc -o progtest $< -lblas -llapack

run: progtest
	./progtest

# Clean up
clean:
	@rm -f progtest