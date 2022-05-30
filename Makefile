# Makefile for ACOTSP-MPI-OMP
# To be adapted to the testbed used

VERSION=MPI-OMP
WITH_ACO=FT_ACO
#HANDLER = FT_ERRORS_ARE_FATAL
HANDLER = FT_ERRORS_RETURN
#HANDLER = FT_ABORT_ON_FAILURE
#HANDLER = FT_IGNORE_ON_FAILURE
OPTIM_FLAGS=-O
WARN_FLAGS=-Wall -Wextra -Werror
# -D_FT_ACO_ to compile with Fault tolerance API
CFLAGS=$(WARN_FLAGS) $(OPTIM_FLAGS) -I ~/mpi402/include/ -fopenmp -DF$(WITH_ACO) -D$(HANDLER)

CC=mpicc
# To change the default timer implementation, uncomment the line below
#TIMER=dos
TIMER=unix
LDLIBS=-lm -L/home/miguel.blanco/mpi402/lib -fopenmp -lmpi -lpthread

all: clean acotsp

clean:
	@$(RM) *.o acotsp

acotsp: acotsp.o parallel.o TSP.o utilities.o ants.o InOut.o $(TIMER)_timer.o ls.o parse.o ft_aco.o

acotsp.o: acotsp.c

TSP.o: TSP.c TSP.h

ants.o: ants.c ants.h

InOut.o: InOut.c InOut.h

utilities.o: utilities.c utilities.h

ls.o: ls.c ls.h

parse.o: parse.c parse.h

parallel.o: parallel.c parallel.h

# Fault tolerance API
ft_aco.o: ft_aco.c ft_aco.h

$(TIMER)_timer.o: $(TIMER)_timer.c timer.h

dist : DIST_SRC_FILES=*.c *.h README *.tsp Makefile gpl.txt
dist : all
	@(mkdir -p ../ACOTSP-$(VERSION)			\
	&& rsync -rlpC --exclude=.svn $(DIST_SRC_FILES) ../ACOTSP-$(VERSION)/ \
        && cd .. 	\
	&& tar cf - ACOTSP-$(VERSION) | gzip -f9 > ACOTSP-$(VERSION).tar.gz \
	&& rm -rf ./ACOTSP-$(VERSION)					\
	&& echo "ACOTSP-$(VERSION).tar.gz created." && cd $(CWD) )
