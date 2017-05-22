/*
 * This library is designed to allow easy implementation of a number
 * of MPI functions. It is primarily based off the book by Peter
 * Pacheco but extensions have been added where they have been found
 * useful.
 */

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <stdbool.h>

#include "MPI_Additions.h"

/****************************************************
 * File Constants                                   *
 ****************************************************/
static int _dest;
static int _source;
static int _my_rank;
static int _my_proc;

static int _my_grid_coord[3]; //Position of _my_rank in grid
static int _grid_dim[3];      //Overal Dimensions of grid
static float _my_local_dim[3];  //Local dimension of the spatial
static float _dim[3];        //Dimensions of system

// region that _my_rank is responsible for

static bool _log;
static char _my_rank_log[100];
static int _my_fd;
static fpos_t _my_pos;

MPI_Status  _status;

/****************************************************
 * MPI Functions                                    *
 ****************************************************/
/* Internal Functions                               */

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

  MPI_Recv(ints          ,
           number        ,
           MPI_INT       ,
           _source       ,
           0             ,
           MPI_COMM_WORLD,
           &_status      );
}

static inline void _sendChar(int number, void ** data){


  printf("rank %d _dest %d number send %d\n",_my_rank,_dest,number);
  fflush( stdout );
  char * chars = (char *)data;

  MPI_Send(chars         ,
           number        ,
           MPI_CHAR      ,
           _dest         ,
           0             ,
           MPI_COMM_WORLD);
}

static inline void _recvChar(int number, void ** data){

  printf("rank %d _source %d number recv %d\n",_my_rank,_source,number);
  fflush( stdout );
  char * chars = (char *) data;

  MPI_Recv(chars         ,
           number        ,
           MPI_CHAR      ,
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
  if(_log){

  }
  printf("_send my_rank %d number %d\n",_my_rank,number);
  fflush( stdout );
  (*fun)(number,data);
  printf("_send my_rank %d number %d\n",_my_rank,number);
  fflush( stdout );
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
  if(_log){

  }
  printf("_recv my_rank %d number %d\n",_my_rank,number);
  fflush( stdout );
  (*fun)(number,data);
  printf("_recv my_rank %d number %d\n",_my_rank,number);
  fflush( stdout );
}

static inline void _send_recv(int number,
                              void ** data,
                              void (*send_fun)(int,void**),
                              void (*recv_fun)(int,void**)){
  if(_my_rank%2==0){
    _send(number,data,send_fun);
  }else{
    _recv(number,data,recv_fun);
  }
}

static inline void _recv_send(int number,
                              void ** data,
                              void (*send_fun)(int,void**),
                              void (*recv_fun)(int,void**)){
  if(_my_rank%2==0){
    _recv(number,data,recv_fun);
  }else{
    _send(number,data,send_fun);
  }
}

static inline void _pingpong(int number,
                             void ** data,
                             void (*send_fun)(int, void **),
                             void (*recv_fun)(int, void **)){

  _dest = _my_rank + 1;
  _source = _my_rank - 1;
  _send_recv(number  ,
             data    ,
             (*send_fun),
             (*recv_fun));

  _source = _my_rank + 1;
  _dest = _my_rank - 1;
  _recv_send(number  ,
             data    ,
             (*send_fun),
             (*recv_fun));
}


/* External Functions                               */

/* test_type 1 - is linear increase in size          *
 * test_type ~1 - is an exponential increase in size */
// res_lin is the resolution 
int pingpong(int test_type   ,
             int max_byte    ,
             int res_lin     ,
             double max_err_perc){

  if(max_byte<0){
    fprintf(stderr,"ERROR max_byte is less than 0\n");
    return -1;
  }

  int power     = 0;
  int num_bytes = 0;
  int inc_bytes = max_byte/res_lin;

  char * bytes = NULL;

  // Linear increase in message size

  while(num_bytes<max_byte){

    double stats[4] = { 0 };
    double err_perc = 100;

    int    iter     = 0;
    // Allocate memory
    bytes = realloc(bytes,sizeof(char)*num_bytes);

    while(err_perc > max_err_perc){
      MPI_Barrier(MPI_COMM_WORLD);

      double t_start = MPI_Wtime();

      printf("rank %d num_bytes %d b\n",_my_rank,num_bytes);
      _pingpong(num_bytes,
          (void **) bytes,
          &_sendChar,
          &_recvChar);

      double t_finish = MPI_Wtime();
      if(_my_rank==0){
        // Divide by 2.0 because there are a total of two messages
        // message one from proc 0 to 1 and then from proc 1 to 0.
        double t_elapsed = (t_finish-t_start)/2.0;
        printf("rank %d num_bytes %d a\n",_my_rank,num_bytes);
        // Finished first iteration of message passing
        iter++;
        calcStats(iter,t_elapsed,stats);
        // Determine if the error in proportion to the mean
        err_perc = stats[3]/stats[1]*100.0;
        MPI_Send(&err_perc, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
      }else{
        MPI_Recv(&err_perc, 1, MPI_DOUBLE, 0, 0,MPI_COMM_WORLD,&_status);
      }
    }

    // Store run in a line in a file
    if(_my_rank==0){
      genStatFile("PerformancePingPong.txt",
                  iter,
                  num_bytes,
                  stats);
    }

    // Linear increase in message size
    if(test_type){
      num_bytes += inc_bytes;
    // else exponential increase in data size
    }else{
      num_bytes = (int) pow(2.0,power);
      power++;
    }

  }

  if(bytes) free(bytes);

  return 0;
}

void init_MPI_Additions(void){
  MPI_Comm_rank(MPI_COMM_WORLD,&_my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&_my_proc);
  sprintf(_my_rank_log,"file_num%d.log",_my_rank);
}

// If given a volume, this function will 
// detemine the optimal grid the processors
// will map to the the volume. GridDim determines what
// dimensiosn we are allowed to map the processor grid too
// if GridDim is 2 even if we have a volume we are mapping 
// to the actual processor grid will be two dimensional 
int calibrateGrid(float Dimensions[3], int GridDim){

  if( Dimensions[0] < 0 || Dimensions[1] < 0 || Dimensions[2] < 0){
    fprintf(stderr,"ERROR Dimension less than 0 in calibrateGrid\n");
    return -1;
  }
  
  {
    int maxGridDim = 3;
    if(Dimensions[0] ==0 ){
      maxGridDim--;
    }
    if(Dimensions[1] == 0){
      maxGridDim--;
    }
    if(Dimensions[2] == 0){
      maxGridDim--;
    }

    if(maxGridDim<GridDim){
      fprintf(stderr,"ERROR maxGridDim allowed is %d you have "
                     "specified a gridDim of %d. Cannot map a "
                     "System of size %g,%g,%g to %d Grid      "
                     "dimensions\n.",maxGridDim,GridDim,
                     Dimensions[0],Dimensions[1],Dimensions[2],
                     GridDim);
      return -1;
    }
  }
  // Determine if any of the dimensions of the physical system are 
  // 0. This will limit the total GridDim allowed. 
  
  // Store the dimensions of the system
  _dim[0] = Dimensions[0];
  _dim[1] = Dimensions[1];
  _dim[2] = Dimensions[2];

  // Determine if lattice can be split
  // up in a 1d 2d or 3d lattice
  int proc_max;
  int proc_mid;
  int proc_min;

  /* Following section determines how best to split up
   * the system, it is assumed that the work will be
   * distributed evenly across the system */
  int * array = malloc(sizeof(int));

  int len_array = primeFactors(_my_proc, &array);

  for(int i=0;i<len_array;i++){
    printf("array[%d] %d\n",i,array[i]);
  }

  int val1 = 1;
  int val2 = 1;
  int val3 = 1;

  /* Determinig how to split up the system based on
   * the system up dimentions and the number
   * of processors. We are trying to keep the size
   * of the blocks as cubic as possible. Cubic volumes
   * have smaller surface are when compared to rectangular
   * cuboid like volumes. */

  /* Here we are making the dimensions as cubic as possible */
  for(int i=len_array-1;i>=0;i--){
    if(val1<=val2 && val1<= val3){
      val1 = array[i]*val1;
    }else if(val2 <= val3 && GridDim>=2){
      val2 = array[i]*val2;
    }else if(GridDim>=3){
      val3 = array[i]*val3;
    }
  }

  int grid_len;
  int grid_wid;
  int grid_hei;

  /* Here we are mapping the processors to the dimensions */
  if(val1>=val2 && val1>= val3){
    proc_max = val1;
    if(val2>=val3){
      proc_mid = val2;
      proc_min = val3;
    }else{
      proc_mid = val3;
      proc_min = val2;
    }
  }else if(val2>=val1 && val2>=val3){
    proc_max = val2;
    if(val1>=val3){
      proc_mid = val1;
      proc_min = val3;
    }else{
      proc_mid = val3;
      proc_min = val1;
    }

  }else{
    proc_max = val3;
    if(val1>= val2){
      proc_mid = val1;
      proc_min = val2;
    }else{
      proc_mid = val2;
      proc_min = val1;
    }
  }

  if(_dim[0]>=_dim[1] && _dim[0]>=_dim[2]){
    grid_len = proc_max;
    if(_dim[1]>=_dim[2]){
      grid_wid = proc_mid;
      grid_hei = proc_min;
    }else{
      grid_hei = proc_mid;
      grid_wid = proc_min;
    }
  }else if(_dim[1]>=_dim[0] && _dim[1]>=_dim[2]){
    grid_wid = proc_max;
    if(_dim[0]>=_dim[2]){
      grid_len = proc_mid;
      grid_hei = proc_min;
    }else{
      grid_hei = proc_mid;
      grid_len = proc_min;
    }
  }else{
    grid_hei = proc_max;
    if(_dim[0]>=_dim[1]){
      grid_len = proc_mid;
      grid_wid = proc_min;
    }else{
      grid_wid = proc_mid;
      grid_len = proc_min;
    }

  }

  _grid_dim[0] = grid_len;
  _grid_dim[1] = grid_wid;
  _grid_dim[2] = grid_hei;

  _my_local_dim[0] = _dim[0] / ((float)_grid_dim[0]);
  _my_local_dim[1] = _dim[1] / ((float)_grid_dim[1]);
  _my_local_dim[2] = _dim[2] / ((float)_grid_dim[2]);

  /* Determine where current rank fits into the grid */
  _my_grid_coord[0] = _my_rank % _grid_dim[0];
  _my_grid_coord[1] = (_my_rank % (_grid_dim[0]*_grid_dim[1]))/_grid_dim[0];
  _my_grid_coord[2] = _my_rank / (_grid_dim[0]*_grid_dim[1]);

  printf("rank %d local_len %f local_wid %f local_hei %f\n",_my_rank,_my_local_dim[0],_my_local_dim[1],_my_local_dim[2]);
  printf("rank %d grid_coord %d %d %d\n",_my_rank,_my_grid_coord[0],_my_grid_coord[1],_my_grid_coord[2]);

  free(array);

  return 0;
}

/****************************************************
 * File IO Functions                                *
 ****************************************************/
/* Internal Functions                               */
/* External Functions                               */

// Function designed to switch from printing to 
// stdout, to a file specific to the rank
void switchStdout(void){
  fflush(stdout);
  fgetpos(stdout,&_my_pos);
  _my_fd = dup(fileno(stdout));
  FILE * fp = freopen(_my_rank_log,"w",stdout);
  if(!fp){
    fprintf(stderr,"Unable to redirect from stdout to file\n");
  }
}

// Function designed to switch from print to 
// a file specific to the rank back to stdout
void revertStdout(void){
  fflush(stdout);
  dup2(_my_fd, fileno(stdout));
  close(_my_fd);
  clearerr(stdout);
  fsetpos(stdout,&_my_pos);
}

/****************************************************
 * Statistics Functions                             *
 ****************************************************/
/* Internal Functions                               */

static inline double _stdDev(double mean, double val, int iter){
  return pow(pow(val-mean,2)/((double)iter+1.0),(1.0/2.0));
}

/* External Functions                               */
int calcStats(int iter,double val,double stats[4]){
  stats[0] += val;                                     // sum
  stats[1]  = stats[0]/(iter+1);                       // running mean
  stats[2]  = _stdDev(stats[1],val,iter);              // standard deviation
  stats[3]  = 1.96*stats[2]/pow((double)iter,1.0/2.0); // error
  return 0;
}

int genStatFile( char * file_name, 
                 int iter, 
                 int num_bytes,
                 double stats[4]){

  FILE *fp;
  fp = fopen(file_name,"a+");
  fprintf(fp,"Iter %8d ",iter);
  fprintf(fp,"Bytes %8d ",num_bytes);
  fprintf(fp,"Mean %12g ",stats[1]);
  fprintf(fp,"Dev %12g ",stats[2]);
  fprintf(fp,"Error+- %12g\n",stats[3]);
  fclose(fp);

  return 0;
}

/****************************************************
 * Math Functions                                   *
 ****************************************************/
/* Internal Functions                               */

/* External Functions                               */

/* Function is designed to print all prime factors  *
 * a given number n                                 */
int primeFactors(int n, int ** array){

  if(!array) return -1;
 
  
    int size = 0;
    if(n==1) {
        (*array)[0] = 1;
        return 1;
    }
    // Print the number of 2s that divide n
    while (n%2 == 0){
        n = n/2;
        if(size>0){
            int * temp = realloc(*array,sizeof(int)*(size+1));
            if(!temp){
                fprintf(stderr,"ERROR temp array is NULL\n");
                exit(1);
            }
            *array = temp;
        }

        (*array)[size] = 2;
        size++;
    }

      // n must be odd at this point.  So we can skip one element
      // (Note i = i +2)
    for (int i = 3; i <= sqrt(n); i = i+2)
    {
        // While i divides n, print i and divide n
        while (n%i == 0)
        {
            n = n/i;
            if(size>0){
                int * temp = realloc(*array,sizeof(int)*(size+1));
                if(!temp){
                    fprintf(stderr,"ERROR temp array is NULL\n");
                    exit(1);
                }
                **array = *temp;
            }
            (*array[size]) = i;
            size++;
        }
    }

    // This condition is to handle the case whien n is a prime number
    // greater than 2
    if (n > 2){
        if(size>0){
            int * temp = realloc(*array,sizeof(int)*(size+1));
            if(!temp){
                fprintf(stderr,"ERROR temp array is NULL\n");
                exit(1);
            }
            **array = *temp;
        }

        (*array)[size] = n;
        size++;
    }
    return size;    // n must be odd at this point.  So we can skip

}



