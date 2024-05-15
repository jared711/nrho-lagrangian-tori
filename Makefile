OPT=-g -Wall
#OPT=-O3 -Wall
CFLAGS=$(OPT) -ffast-math -fdiagnostics-color=always

UTILS=param

all : $(UTILS)

###########
# Utilites
###########s
# - Secció de Poincaré amb RTBP
param : param.cc seccp.o fluxvp.o rk78vp.o rtbphp.o campvp.o scread.o vbprintf.o
	g++ -o param $(OPT) $(CFLAGS) param.cc seccp.o fluxvp.o rk78vp.o rtbphp.o campvp.o scread.o vbprintf.o -lm

##########
# Objects
##########
seccp.o : seccp.c
	gcc -c $(CFLAGS) $<

rtbphp.o : rtbphp.c
	gcc -c $(CFLAGS) $<

fluxvp.o : fluxvp.c
	gcc -c $(CFLAGS) $<

scread.o : scread.c
	gcc -c $(CFLAGS) $<

vbprintf.o : vbprintf.c
	gcc -c $(CFLAGS) $<

rk78vp.o : rk78vp.c
	gcc -c $(CFLAGS) $<

campvp.o : campvp.c
	gcc -c $(CFLAGS) $<

# matrix.o : matrix.h
# 	gcc -c $(CFLAGS) $<

# grid.o : grid.h
# 	gcc -c $(CFLAGS) $<

# complex.o : complex.h
# 	gcc -c $(CFLAGS) $<

###############
# Housekeeping
###############
clean :
	rm -f *.o
realclean : clean
	rm -f $(UTILS)
