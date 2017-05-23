

/* This function is called after MPI_Init is used *
 * it updates variables with in the library that  *
 * have file scope.                               */
void init_MPI_Additions(void);

int getMyRank(void);

int getMyProc(void);

char * getMyRankLog(void);

/* The pingpong function is used to determine the *
 * message passing performance of a parallel      *
 * system by calling MPI_Send and MPI_Recv. The   *
 * variables are defined as:                      *
 * test_type - This can be either 0 or a non zero *
 *             number. 0 indicates a linearly     *
 *             increasing message is used. A non  *
 *             zero number indicates that an exp  *
 *             increasing message size is used.   *
 * max_byte  - Determines the maximum size of the *
 *             message to be passed.              *
 * res_lin   - If using a linear test determines  *
 *             the resolution of the messages that*
 *             are passed. I.e. The increment size*
 *             will be defined by max_byte/res_lin*
 * max_err_perc - Defines the maximum amount of   *
 *                error accepted. Will do multiple*
 *                iterations of the pingpong until*
 *                the error is less than the max  */
int pingpong(int test_type,
             int max_byte ,
             int res_lin  ,
             double max_err_perc);



/* This function works by determining all of the  *
 * prime factors that exist within the interger n.*
 * The number of prime numbers is then returned in*
 * the array.                                     */
int primeFactors(int n, int ** array);

/* This function works by updating the statistics *
 * stores in the stats array.                     */
int calcStats(int iter,
              double val,
              double stats[4]);

/* This function will right a line of data to the *
 * file specified by file_name. The statistics    *
 * included are:                                  *
 * iter      - the iteration number specified for *
 *             each data point for instance if the*
 *             error is to large with one data    *
 *             point we will use two thus iter    *
 *             will have a value of two etc...    *
 * num_bytes - number of bytes                    *
 * stats[1] - the running mean                    *
 * stats[2] - the standard deviation              *
 * stats[3] - the Error                           *
 *                                                *
 * Note that stat[0] the sum is not used          */
int genStatFile( char * file_name,
                 int iter,
                 int num_bytes,
                 double stats[4]);

/* This function works by mapping a physical     *
 * to the number of processors. The idea behind  *
 * the algorithm is to achieve optimal           *
 * performance we want to maximize the volume to *
 * surface area. Hence the physical volume is    *
 * divided up into rectoids that as closely      *
 * resemble cubes as possible. The user is       *
 * however, allowed to specify the dimensions    *
 * the physical volume can be divided into. So if*
 * if I have a physical volume of 3 dimensions I *
 * can specify that I want it mapped to a        *
 * processor grid of 2 dimensions. The parameters*
 * are defined as follows:                       *
 * Dimensions - The physical dimensiosn of the   *
 *              system.                          *
 * GridDim    - The dimensions of the processor  *
 *              grid that the pysical volume will*
 *              be mapped to.                    *
 *                                               *
 * NOTE: If the physical system is only 2d the   *
 * processor grid can not have 3 dimensions      *
 * NOTE: GridDim is only allowed to have values  *
 * between 1 and 3                               */
int calibrateGrid( float Dimensions[3], int GridDim);

/* This function is useful when debugging using  *
 * printf functions it will switch from printing *
 * to stdout to a log file specific to each      *
 * processor.                                    */
void switchStdout(void);

/* To revert the behavior and have stdout print  *
 * to the terminal again the following function  *
 * can be called.                                */
void revertStdout(void);

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
                              int const num_int);

/* This function works once the processor grid *
 * has been calibrated. For a set of           *
 * within the processor grid it will get the   *
 * appropriate rank.                           */
int get_proc_rank(int grid_coord_x,int grid_coord_y, int grid_coord_z);
