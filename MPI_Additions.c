/*
 * This library is designed to allow easy implementation of a number
 * of MPI functions. It is primarily based off the book by Peter
 * Pacheco but extensions have been added where they have been found
 * useful.
 */

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#inlcude "MPI_Additions"

static int _dest;
static int _source;
static int _my_rank;
static int _my_proc;

static bool _log;
static char _my_rank_log[100];

MPI_Status  _status;

int init_MPI_Additions(void){
  MPI_Comm_rank(MPI_COMM_WORLD,&_my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&_my_proc);
}

static inline double _stdDev(double mean, int iter){

}

static inline void _sendDouble(int number, void ** data){

  double * doubles = (double *)data;

  MPI_Send(doubles       ,
           number        ,
           MPI_DOUBLE    ,
           _dest         ,
           0             ,
           MPI_COMM_WORLD);
}

static inline void _recvDouble(int number, void ** data){

  double * doubles = (double *) data;

  MPI_Recv(doubles       ,
           number        ,
           MPI_DOUBLE    ,
           _source       ,
           0             ,
           MPI_COMM_WORLD,
           &_status      );
}

static inline void _sendInt(int number, void ** data){

  int * ints = (int *)data;

  MPI_Send(ints          ,
           number        ,
           MPI_INT       ,
           _dest         ,
           0             ,
           MPI_COMM_WORLD);
}

static inline void _recvInt(int number, void ** data){

  int * ints = (int *) data;

  MPI_Recv(doubles       ,
           number        ,
           MPI_INT    ,
           _source       ,
           0             ,
           MPI_COMM_WORLD,
           &_status      );
}

static inline void _sendChar(int number, void ** data){

  char * chars = (char *)data;

  MPI_Send(chars         ,
           number        ,
           MPI_CHAR      ,
           _dest         ,
           0             ,
           MPI_COMM_WORLD);
}

static inline void _recvChar(int number, void ** data){

  int * ints = (int *) data;

  MPI_Recv(doubles       ,
           number        ,
           MPI_INT    ,
           _source       ,
           0             ,
           MPI_COMM_WORLD,
           &_status      );
}


static inline void _send(int number, void ** data, void (*fun)(int, void **)){
  #ifdef _ERROR_CHECKING_ON_
  if(number<0){
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return;
  }
  #endif
  if(log){

  }

  *fun(number,data);
}

static inline void _recv(int number, void ** data, void (*fun)(int, void **)){
  #ifdef _ERROR_CHECKING_ON_
  if(number<0){
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return;
  }
  #endif
  if(log){

  }
  *fun(number,data);
}

static inline void _send_recv(int number,
                              void ** data,
                              void (*send_fun)(int,void**),
                              void (*recv_fun)(int,void**)){
  if(my_rank%2==0){
    _send(number,data,send_fun);
  }else{
    _recv(number,data,recv_fun);
  }
}

static inline void _recv_send(int number,
                              void ** data,
                              void (*send_fun)(int,void**),
                              void (*recv_fun)(int,void**)){
  if(my_rank%2==0){
    _recv(number,data,recv_fun);
  }else{
    _send(number,data,send_fun);
  }
}

static inline void _pingpong(int number,
                             void ** data,
                             void (*send_fun)(int, void **)
                             void (*recv_run)(int, void **)){

  _dest = _my_rank + 1;
  _source = _my_rank - 1;
  _send_recv(number  ,
             data    ,
             send_fun,
             recv_fun);

  _source = _my_rank+1;
  _dest = _my_rank-1;
  _recv_send(number  ,
             data    ,
             send_fun,
             recv_fun);
}

int pingpong(int test_type,
             int max_byte ,
             int res_lin  ){

  if(max_byte<0){
    fprintf(stderr,"ERROR max_byte is less than 0\n");
    return -1;
  }

  int power     = 0;
  int num_bytes = 0;
  int inc_bytes = max_byte/res_lin;

  // Linear increase in message size

  while(num_bytes<max_byte){

    // Allocate memory
    char * bytes = calloc(bytes,sizeof(char));

    MPI_Barrier(MPI_COMM_WORLD);

    double t_start = MPI_Wtime();

    _pingpong(number_bytes,
              data,
              &sendChar,
              &recvChar);

    double t_finish = MPI_Wtime();

    double t_elapsed = (finish-start)/2.0;

    // Linear increase in message size
    if(test_type){
      num_bytes += inc_bytes;
    // else exponential increase in data size
    }else{
      num_bytes = (int) pow(2.0,power)
      power++;
    }
  }
  return 0;
}
