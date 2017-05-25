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

static int   _my_grid_coord[3]; //Position of _my_rank in grid
static int   _grid_dim[3];      //Overal Dimensions of grid
static float _my_local_dim[3];  //Local dimension of the spatial
static float _dim[3];           //Dimensions of system

// region that _my_rank is responsible for

static bool   _log;
static char   _my_rank_log[100];
static int    _my_fd;
static fpos_t _my_pos;

MPI_Status  _status;
MPI_Group   _world_group;

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

  char * chars = (char *)data;

  MPI_Send(chars         ,
           number        ,
           MPI_CHAR      ,
           _dest         ,
           0             ,
           MPI_COMM_WORLD);
}

static inline void _recvChar(int number, void ** data){

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
  (*fun)(number,data);
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
  (*fun)(number,data);
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

  #ifdef _ERROR_CHECKING_ON_
  if(_my_proc<2){
    fprintf(stderr,"ERROR need at least two processors to"
                   " use the pingpong function\n");
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(max_byte<0){
    fprintf(stderr,"ERROR max_byte is less than 0\n");
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(res_lin<2 && test_type){
    fprintf(stderr,"ERROR linear resolution or res_lin must"
                   " be greater than two when test_type is "
                   "not 0.\n");
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(max_err_perc<0){
    fprintf(stderr,"ERROR max error percent must be a non "
                   "negative number\n");
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  // Create a Group consisting of only the first two processors
  const int ranks[2] = {0, 1};
  MPI_Group group_01;
  MPI_Group_incl(_world_group,2,ranks,&group_01);
  // Create a new communicator based on the group
  MPI_Comm MPI_01_COMM;
  MPI_Comm_create_group(MPI_COMM_WORLD,group_01,0,&MPI_01_COMM);


  // Only runs with the first two processors
  if(_my_rank<2){
    
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
        // Only block the two processors in the 01 Communicator
        MPI_Barrier(MPI_01_COMM);

        double t_start = MPI_Wtime();

        _pingpong(num_bytes,
            (void **) bytes,
            &_sendChar,
            &_recvChar);

        double t_finish = MPI_Wtime();
        if(_my_rank==0){
          // Divide by 2.0 because there are a total of two messages
          // message one from proc 0 to 1 and then from proc 1 to 0.
          double t_elapsed = (t_finish-t_start)/2.0;
          // Finished first iteration of message passing
          iter++;
          calcStats(iter,t_elapsed,stats);
          // Determine if the error in proportion to the mean
          err_perc = stats[3]/stats[1]*100.0;
          MPI_Send(&err_perc, 1, MPI_DOUBLE, 1, 0, MPI_01_COMM);
        }else{
          MPI_Recv(&err_perc, 1, MPI_DOUBLE, 0, 0,MPI_01_COMM,&_status);
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
  }

  MPI_Group_free(&group_01);
  MPI_Comm_free(&MPI_01_COMM);

  return 0;
}

void init_MPI_Additions(void){
  MPI_Comm_rank(MPI_COMM_WORLD,&_my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&_my_proc);
  MPI_Comm_group(MPI_COMM_WORLD,&_world_group);
  sprintf(_my_rank_log,"file_num%d.log",_my_rank);
}

void clean_MPI_Additions(void){
  MPI_Group_free(&_world_group);
}

int getMyRank(void){
  return _my_rank;
}

int getMyGridX(void){
  return _my_grid_coord[0];
}

int getMyGridY(void){
  return _my_grid_coord[1];
}

int getMyGridZ(void){
  return _my_grid_coord[2];
}

int getMyProc(void){
  return _my_proc;
}

char * getMyRankLog(void){
  return &_my_rank_log[0];
}

// If given a volume, this function will 
// detemine the optimal grid the processors
// will map to the the volume. GridDim determines what
// dimensiosn we are allowed to map the processor grid too
// if GridDim is 2 even if we have a volume we are mapping 
// to the actual processor grid will be two dimensional 
int calibrateGrid(float Dimensions[3], int GridDim){

  #ifdef _ERROR_CHECKING_ON_
  if( Dimensions[0] < 0 || Dimensions[1] < 0 || Dimensions[2] < 0){
    fprintf(stderr,"ERROR Dimension less than 0 in calibrateGrid\n");
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(GridDim<1){
    fprintf(stderr,"ERROR lowest allowed grid dimension is 1\n");
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
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
      #ifdef _FORCE_HARD_CRASH_
      exit(1);
      #endif
      return -1;
    }
  }
  #endif
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
 * a given number n, and return them to the int **  *
 * ptr designated by array. This function will      *
 * return the size of the array upon completion.    *
 * array must be created using dynamically allocated*
 * memory for the size of one int i.e.              *
 *                                                  *
 * int * array = malloc(sizeof(int));               *
 * primeFactors(5,&array);                          *
 *                                                  */
int primeFactors(int n, int ** array){

  #ifdef _ERROR_CHECKING_ON_
  if(!array) {
    fprintf(stderr,"ERROR array is NULL\n");
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(n<1){
    fprintf(stderr,"ERROR n must be a non negative "
                   "integer\n");
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  int size = 0;
  int capacity = 1;
  (*array)[0] = 1;
  size++;

  if(n==1) {
    return 1;
  }
  // Print the number of 2s that divide n
  while (n%2 == 0){
    n = n/2;
    if(size>=capacity){
      capacity = size*2+1;
      int * temp = realloc(*array,sizeof(int)*capacity);
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
      if(size>=capacity){
        capacity = size*2+1;
        int * temp = realloc(*array,sizeof(int)*capacity);
        if(!temp){
          fprintf(stderr,"ERROR temp array is NULL\n");
          exit(1);
        }
        **array = *temp;
      }
      (*array)[size] = i;
      size++;
    }
  }

  // This condition is to handle the case whien n is a prime number
  // greater than 2
  if (n > 2){
    if(size>=capacity){
      capacity = size*2+1;
      int * temp = realloc(*array,sizeof(int)*capacity);
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


/* Passes integers between neighbors
 * hig is passed to the neighbor above
 * low is passed to the neighbor below
 * the value of hig is received by the
 * integer neigh_hig
 * the value of low is received by the
 * integer neigh_low
 * Direction specfies the direction
 * that the communication will take
 * place where:
 * 0 - x direction
 * 1 - y direction
 * 2 - z direction
 * num_int - number of integers */
int MPI_Communicate_Neighbors(int const * low,
                              int const * hig,
                              int * neigh_low,
                              int * neigh_hig,
                              int const Direction,
                              MPI_Request *request_l_recv,
                              MPI_Request *request_l_send,
                              MPI_Request *request_h_recv,
                              MPI_Request *request_h_send,
                              int const num_int){

  if(Direction>2 || Direction<0){
    fprintf(stderr,"ERROR Direction must be a number between 0 and 2:\n"
                       "0 - x dim\n"
                       "1 - y dim\n"
                       "2 - z dim\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  /************************************/
  /* Even processors communicate with
   * odd processors to the right */

  /* 0 -> 1   2 ->
   * Processor 0 must receive from 2 */
  int neigh;
  int tag = 0;

  int change_x = 0;
  int change_y = 0;
  int change_z = 0;
  if(Direction==0){
    change_x++;
  }else if(Direction==1){
    change_y++;
  }else if(Direction==2){
    change_z++;
  }

  fflush(stdout);
  printf("rank %d up_x %d up_y %d up_z %d\n",_my_rank,
      get_proc_rank(_my_grid_coord[0]+1,_my_grid_coord[1],_my_grid_coord[2])==MPI_PROC_NULL ? _my_rank : get_proc_rank(_my_grid_coord[0]+1,_my_grid_coord[1],_my_grid_coord[2]),
      get_proc_rank(_my_grid_coord[0],_my_grid_coord[1]+1,_my_grid_coord[2])==MPI_PROC_NULL ? _my_rank : get_proc_rank(_my_grid_coord[0],_my_grid_coord[1]+1,_my_grid_coord[2]),
      get_proc_rank(_my_grid_coord[0],_my_grid_coord[1],_my_grid_coord[2]+1)==MPI_PROC_NULL ? _my_rank : get_proc_rank(_my_grid_coord[0],_my_grid_coord[1],_my_grid_coord[2]+1));
  fflush(stdout);
  printf("_grid_dim %d %d %d\n",_grid_dim[0],_grid_dim[1],_grid_dim[2]);
  if(_my_grid_coord[Direction]%2==0){
    /* Even Send to right */
    neigh = get_proc_rank(_my_grid_coord[0]+change_x,
        _my_grid_coord[1]+change_y,
        _my_grid_coord[2]+change_z);

    printf("my_rank %d send %d tag %d hig %d\n",_my_rank,
neigh,tag,*hig);
    //MPI_Send(&hig,1,MPI_INT,neigh,tag,MPI_COMM_WORLD);
    MPI_Isend(hig,num_int,MPI_INT,
        neigh,
        tag,
        MPI_COMM_WORLD,
        request_h_send);
    /* If odd number of dimensions in the x grid direction
     * rank 0 must also receive after sending */
    if(_grid_dim[Direction]%2!=0 && _my_grid_coord[Direction]==0){
      /* 0 receive from left */
      neigh = get_proc_rank(_my_grid_coord[0]-change_x,
          _my_grid_coord[1]-change_y,
          _my_grid_coord[2]-change_z);

      printf("my_rank %d recv %d tag %d\n",_my_rank, neigh,tag);
      MPI_Irecv(neigh_hig,num_int,MPI_INT,
          neigh,
          tag,
          MPI_COMM_WORLD,
          request_h_recv);
      //MPI_Recv(&neigh_hig,1,MPI_INT,neigh,tag,MPI_COMM_WORLD,&status);
      printf("my_rank %d receiving from neigh %d neigh_hig, %d\n",_my_rank,neigh,*neigh_hig);
    }
  }else{
    /* Odd receive from left */
    neigh = get_proc_rank(_my_grid_coord[0]-change_x,
        _my_grid_coord[1]-change_y,
        _my_grid_coord[2]-change_z);

    printf("my_rank %d recv %d tag %d\n",_my_rank, neigh,tag);
    MPI_Irecv(neigh_hig,num_int,MPI_INT,
        neigh,
        tag,
        MPI_COMM_WORLD,
        request_h_recv);
    //MPI_Recv(&neigh_hig,1,MPI_INT,neigh,tag,MPI_COMM_WORLD,&status);
    printf("my_rank %d receiving from neigh %d neigh_hig, %d\n",_my_rank,neigh, *neigh_hig);
  }

  /* Send to the other neighbor on the x side */
  /* 0 <- 1   2 <-
   * Processor 0 must receive from 2 */
  tag = 1;
  if(_my_grid_coord[Direction]%2==0){
    /* Even receive from right */
    neigh = get_proc_rank(_my_grid_coord[0]+change_x,
        _my_grid_coord[1]+change_y,
        _my_grid_coord[2]+change_z);
    printf("my_rank %d recv %d tag %d\n",_my_rank, neigh,tag);
    MPI_Irecv(neigh_low,num_int,MPI_INT,
        neigh,
        tag,
        MPI_COMM_WORLD,
        request_l_recv);
    //MPI_Recv(&neigh_low,1,MPI_INT,neigh,tag,MPI_COMM_WORLD,&status);
    printf("my_rank %d receiving from neigh %d neigh_low, %d\n",_my_rank,neigh,*neigh_low);

    /* If odd number of dimensions in the x grid direction
     * rank 0 must also send after receiving */
    if(_grid_dim[Direction]%2!=0 && _my_grid_coord[Direction]==0){
      /* 0 send to left*/
      neigh = get_proc_rank(_my_grid_coord[0]-change_x,
          _my_grid_coord[1]-change_y,
          _my_grid_coord[2]-change_z);
      printf("my_rank %d send %d tag %d low %d\n",_my_rank,
neigh,tag,*low);
      MPI_Isend(low,num_int,MPI_INT,
          neigh,
          tag,
          MPI_COMM_WORLD,
          request_l_send);
      //MPI_Send(&low,1,MPI_INT,neigh,tag,MPI_COMM_WORLD);
      printf("my_rank %d sending low %d\n",_my_rank, *low);

    }
  }else{
    /* Odd send to left */
    neigh = get_proc_rank(_my_grid_coord[0]-change_x,
        _my_grid_coord[1]-change_y,
        _my_grid_coord[2]-change_z);

    printf("my_rank %d send %d tag %d low %d\n",_my_rank,
neigh,tag,*low);
    MPI_Isend(low,num_int,MPI_INT,
        neigh,
        tag,
        MPI_COMM_WORLD,
        request_l_send);
    //MPI_Send(&low,1,MPI_INT,neigh,tag,MPI_COMM_WORLD);
    printf("my_rank %d sending low %d\n",_my_rank,*low);
  }

  /************************************/
  /* Even processors communicate with
   * odd processors to the left */
  /* Now we must communicate in the other direction */
  fflush(stdout);
  /* 0    1 <- 2 */
  /* If even
   * 0    1 <- 2    3 <-*/
  tag = 4;
  if(_grid_dim[Direction]%2==1){
    if(_my_grid_coord[Direction]%2==0){
      if(_my_grid_coord[Direction]!=0){
        /* Even Send to left with the exception of proc-
         * essor 0 which has alread sent */
        neigh = get_proc_rank(_my_grid_coord[0]-change_x,
            _my_grid_coord[1]-change_y,
            _my_grid_coord[2]-change_z);
        printf("my_rank %d send %d \n",_my_rank, neigh);
        fflush(stdout);
        MPI_Isend(low,num_int,MPI_INT,
            neigh,
            tag,
            MPI_COMM_WORLD,
            request_l_send);
        //MPI_Send(&low,1,MPI_INT,neigh,tag,MPI_COMM_WORLD);
        printf("my_rank %d sending low %d\n",_my_rank,*low);
      }
    }else{
      /* Odd receive from the right */
      neigh = get_proc_rank(_my_grid_coord[0]+change_x,
          _my_grid_coord[1]+change_y,
          _my_grid_coord[2]+change_z);

      printf("my_rank %d recv %d \n",_my_rank, neigh);
      fflush(stdout);
      MPI_Irecv(neigh_low,num_int,MPI_INT,
          neigh,
          tag,
          MPI_COMM_WORLD,
          request_l_recv);
      //MPI_Recv(&neigh_low,1,MPI_INT,neigh,tag,MPI_COMM_WORLD,&status);
      printf("my_rank %d recv from neigh %d neigh_low %d\n",_my_rank,neigh,*neigh_low);
    }
  }else{
    if(_my_grid_coord[Direction]%2==0){
      /* Even Send to left with the exception of proc-
       * essor 0 which has alread sent */
      neigh = get_proc_rank(_my_grid_coord[0]-change_x,
          _my_grid_coord[1]-change_y,
          _my_grid_coord[2]-change_z);
      printf("my_rank %d send %d \n",_my_rank, neigh);
      fflush(stdout);
      MPI_Isend(low,num_int,MPI_INT,
          neigh,
          tag,
          MPI_COMM_WORLD,
          request_l_send);
      //MPI_Send(&low,1,MPI_INT,neigh,tag,MPI_COMM_WORLD);
      printf("my_rank %d sending low %d\n",_my_rank,*low);
    }else{
      /* Odd receive from the right */
      neigh = get_proc_rank(_my_grid_coord[0]+change_x,
          _my_grid_coord[1]+change_y,
          _my_grid_coord[2]+change_z);

      printf("my_rank %d recv %d \n",_my_rank, neigh);
      fflush(stdout);
      MPI_Irecv(neigh_low,num_int,MPI_INT,
          neigh,
          tag,
          MPI_COMM_WORLD,
          request_l_recv);
      //MPI_Recv(&neigh_low,1,MPI_INT,neigh,tag,MPI_COMM_WORLD,&status);
      printf("my_rank %d recv from neigh %d neigh_low %d\n",_my_rank,neigh,*neigh_low);
    }


  }
  /* Odd number of grid dimensions */
  /* 0    1 -> 2 */
  /* Even number of grid dimensions
   * 0    1 -> 2   3 -> */
  tag = 5;
  if(_grid_dim[Direction]%2==1){
    if(_my_grid_coord[Direction]%2!=0){
      /* Odd send to right */
      neigh = get_proc_rank(_my_grid_coord[0]+change_x,
          _my_grid_coord[1]+change_y,
          _my_grid_coord[2]+change_z);

      printf("my_rank %d send %d \n",_my_rank, neigh);
      fflush(stdout);
      MPI_Isend(hig,num_int,MPI_INT,
          neigh,
          tag,
          MPI_COMM_WORLD,
          request_h_send);
      //MPI_Send(&hig,1,MPI_INT,neigh,tag,MPI_COMM_WORLD);
    }else{

      if(_my_grid_coord[Direction]!=0){
        /* Even receive from the left with the exception
         * of processor 0 which has already received from
         * the left if odd number of processors */
        neigh = get_proc_rank(_my_grid_coord[0]-change_x,
            _my_grid_coord[1]-change_y,
            _my_grid_coord[2]-change_z);

        printf("my_rank %d recv %d \n",_my_rank, neigh);
        fflush(stdout);
        MPI_Irecv(neigh_hig,num_int,MPI_INT,
            neigh,
            tag,
            MPI_COMM_WORLD,
            request_h_recv);
        //MPI_Recv(&neigh_hig,1,MPI_INT,neigh,tag,MPI_COMM_WORLD,&status);
      }
    }
  }else{
    if(_my_grid_coord[Direction]%2!=0){
      /* Odd send to right */
      neigh = get_proc_rank(_my_grid_coord[0]+change_x,
          _my_grid_coord[1]+change_y,
          _my_grid_coord[2]+change_z);

      printf("my_rank %d send %d \n",_my_rank, neigh);
      fflush(stdout);
      MPI_Isend(hig,num_int,MPI_INT,
          neigh,
          tag,
          MPI_COMM_WORLD,
          request_h_send);
      //MPI_Send(&hig,1,MPI_INT,neigh,tag,MPI_COMM_WORLD);
    }else{

      /* Even receive from the left with the exception
       * of processor 0 which has already received from
       * the left if odd number of processors */
      neigh = get_proc_rank(_my_grid_coord[0]-change_x,
          _my_grid_coord[1]-change_y,
          _my_grid_coord[2]-change_z);

      printf("my_rank %d recv %d \n",_my_rank, neigh);
      fflush(stdout);
      MPI_Irecv(neigh_hig,num_int,MPI_INT,
          neigh,
          tag,
          MPI_COMM_WORLD,
          request_h_recv);
      //MPI_Recv(&neigh_hig,1,MPI_INT,neigh,tag,MPI_COMM_WORLD,&status);
    }

  }

  printf("neigh_hig %d neigh_low %d\n",*neigh_hig,*neigh_low);
  return 0;
}

int get_proc_rank(int grid_coord_x,int grid_coord_y, int grid_coord_z){
  /* Apply Periodic conditions */
  grid_coord_x = (grid_coord_x+_grid_dim[0])%_grid_dim[0];
  grid_coord_y = (grid_coord_y+_grid_dim[1])%_grid_dim[1];
  grid_coord_z = (grid_coord_z+_grid_dim[2])%_grid_dim[2];

  int val = _grid_dim[0]*_grid_dim[1]*grid_coord_z+_grid_dim[0]*grid_coord_y+grid_coord_x;
  printf("grid_coord_x %d grid_coord_y %d grid_coord_z %d dimx %d dimy %d dimz %d rank %d\n",
      grid_coord_x, grid_coord_y,grid_coord_z,_grid_dim[0],_grid_dim[1],_grid_dim[2],val);

  if(val==_my_rank){
    val = MPI_PROC_NULL;
  }

  return val;
}
