#include <iostream>
#include "mpi.h"
#include <vector>
int main(int argc, char* argv[]) {

  double* localB; // Те составляющие, которые необходимы для вычисления(приходят 1 раз)
  int rank;
  int countProc;
  //Блок инициализации
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &countProc);


  MPI_Finalize();
}
