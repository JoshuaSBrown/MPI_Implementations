#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "MPI_Additions.h"

int main(void){

  MPI_Init(NULL, NULL);
  printf("Testing: init_MPI_Additions\n");
  init_MPI_Additions();
  printf("Rank: %d\n",getMyRank());
  printf("Proc: %d\n",getMyProc());
  printf("Log:  %s\n",getMyRankLog());

  printf("Testing: pingpong\n");
  if(getMyRank()<2){
    int rv = pingpong(1,10000,20,5.0);
    assert(rv==0);
    rv = pingpong(1,-1,20,5.0);
    assert(rv==-1);
    rv = pingpong(1,10000,1,5.0);
    assert(rv==-1);
    rv = pingpong(1,10000,1,-1.0);
    assert(rv==-1);
  }

  printf("Testing: switchStdout & revertStdout\n");
  {
    switchStdout();
    printf("Testing 1\n");
    revertStdout();
    printf("Testing 2\n"); 
  }

  printf("Testing: primeFactors\n");
  {
    int * array = malloc(sizeof(int));
    int rv = primeFactors(0,&array);
    assert(rv==-1);
    rv = primeFactors(1,NULL);
    assert(rv==-1);
    rv = primeFactors(1,&array);
    assert(rv==1);
    assert(array[0]==1);
    rv = primeFactors(2,&array);
    assert(rv==2);
    assert(array[0]==1);
    assert(array[1]==2);
    rv = primeFactors(3,&array);
    assert(rv==2);
    assert(array[0]==1);
    assert(array[1]==3);
    rv = primeFactors(4,&array);
    assert(rv==3);
    assert(array[0]==1);
    assert(array[1]==2);
    assert(array[2]==2);
    rv = primeFactors(5,&array);
    assert(rv==2);
    assert(array[0]==1);
    assert(array[1]==5);
    rv = primeFactors(6,&array);
    assert(rv==3);
    assert(array[0]==1);
    assert(array[1]==2);
    assert(array[2]==3);
    free(array);
  }

  printf("Testing: calibrateGrid\n");
  if(getMyProc()!=4){
    printf("Will not run calibrateGrid test because you must specify 4 processors and you have specified %d processors\n",getMyProc());
  }else{
    float Dims[3] = { 100.0, 100.0, 100};
    calibrateGrid( Dims, 2);
    if(getMyRank()==0){
      int x = getMyGridX();
      int y = getMyGridY();
      int z = getMyGridZ();
      assert(x==0);
      assert(y==0);
      assert(z==0);
    }else if(getMyRank()==1){
      int x = getMyGridX();
      int y = getMyGridY();
      int z = getMyGridZ();
      assert(x==1);
      assert(y==0);
      assert(z==0);
    }else if(getMyRank()==2){
      int x = getMyGridX();
      int y = getMyGridY();
      int z = getMyGridZ();
      assert(x==0);
      assert(y==1);
      assert(z==0);
    }else if(getMyRank()==3){
      int x = getMyGridX();
      int y = getMyGridY();
      int z = getMyGridZ();
      assert(x==1);
      assert(y==1);
      assert(z==0);
    }
  }

  MPI_Finalize();
  return 0;
}
