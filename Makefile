CC = mpicc
CFLAGS = -Wall -Wextra -std=c99

ALL : test_MPI_Additions

test_MPI_Additions : test_MPI_Additions.c MPI_Additions.o
	$(CC) $(CFLAGS) -o test_MPI_Additions MPI_Additions.o

MPI_Additions.o : MPI_Additions.c
	$(CC) $(CFLAGS) -c MPI_Additions.c

.PHONY : clean
clean :
		$(RM) *.o test_MPI_Additions
