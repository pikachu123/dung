CC		= g++
hBOALIBDIR	= ..

CPPFLAGS	= -I$(hBOALIBDIR) -c -Wall  -g#-O3
LFLAGS		=

all:	run_hBOA

clean:
	-rm -f *.o
	-rm -f *.out.fitness.*
	-rm -f run_hBOA
	make -C $(hBOALIBDIR) clean

run_hBOA: run_hBOA.o $(hBOALIBDIR)/libhBOA.a
	$(CC) $(LFLAGS) -o $@ $^ 

$(hBOALIBDIR)/libhBOA.a:
	make -C $(hBOALIBDIR)
