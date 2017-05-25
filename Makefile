CC = mpicc
CFLAGS = -Wall -Wextra -std=gnu99

ALL : test_MPI_Additions

test_MPI_Additions : test_MPI_Additions.c MPI_Additions.o
	$(CC) $(CFLAGS) -o test_MPI_Additions test_MPI_Additions.c MPI_Additions.o -lm

MPI_Additions.o : MPI_Additions.c
	$(CC) $(CFLAGS) -c MPI_Additions.c -lm

.PHONY : clean
clean :
		$(RM) *.o test_MPI_Additions PerformancePingPong.txt file_num*.log
