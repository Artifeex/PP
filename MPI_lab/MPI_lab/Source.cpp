#include <iostream>
#include "mpi.h"
#include <vector>
int main(int argc, char* argv[]) {

  double* localB; // �� ������������, ������� ���������� ��� ����������(�������� 1 ���)
  int rank;
  int countProc;
  //���� �������������
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &countProc);


  MPI_Finalize();
}
