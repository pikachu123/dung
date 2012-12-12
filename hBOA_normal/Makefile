#  =======================================================
#
#  Makefile for the Bayesian optimization Algorithm (BOA)
#
#  author:        Martin Pelikan
#  last modified: February, 1999
#
#  modified by Ying-ping Chen
#  November 13, 2002
#  for Linux binary library distribution
#
#  =======================================================

AR		= ar
RANLIB		= ranlib

CPPFLAGS	= -c -Wall  -g  -rdynamic#-O3
LFLAG		=

OBJS  = args.o		\
	arrays.o	\
	BDe.o           \
	bic.o           \
	binary.o	\
	boa.o           \
	copy.o		\
        counts.o        \
	decisionGraph.o  \
	discretization.o      \
	decisionGraphBoa.o\
        distance.o      \
	ea.o		\
	fitness.o	\
	frequencyDecisionGraph.o\
	frequencyTree.o \
	getFileArgs.o	\
	graph.o		\
	hBoa.o          \
	hBOAmain.o	\
	index.o		\
	individual.o	\
	labeledTreeNode.o\
        mdl.o           \
	memalloc.o	\
	mymath.o        \
	population.o	\
	priors.o	\
	random.o	\
	recombination.o	\
	reordering.o	\
	replace.o	\
	select.o	\
	stack.o		\
	statistics.o	\
	userBOA.o	\
	utils.o		\

#
# Create libhBOA.a
#
libhBOA.a: $(OBJS)
	$(AR) -rc $@ $^
	$(RANLIB) $@

#
# make clean removes all .o files
#

clean: 
	@rm -f $(OBJS)
	@rm -f libhBOA.a

# make all recompiles everything 

all:
	@make clean
	@make libhBOA.a
