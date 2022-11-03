### snapshot options #######
EXTRAS +=

#CC
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=4
DC     += -DLOCK
GSLL   := -lgsl -lgslcblas
CFLAGS := -Wall -O3 -march=native -ftree-vectorize -fopenmp -g
LIBS   := -lm $(GSLL) 

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS :=  variables.o io.o allocate.o leesnap.o grid.o iden.o 

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := main.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

main.x: main.c $(OBJS)
	$(CC) $(CFLAGS) $(EXTRAS) $^  -o $@ $(LIBS)

clean:
	rm -rf $(OBJS)
	rm -rf main.o
	rm -rf $(EXEC)
