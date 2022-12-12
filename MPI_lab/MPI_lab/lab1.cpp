#include "mpi.h"
#include <iostream>
#include <random>
#include <ctime>

int main(int argc, char* argv[])
{
  int rank;
  int remainder;
  int countProc;
  int* matrix = nullptr;
  int countRow = 10;
  int lenRow;
  int* result = nullptr;
  int* localResult = nullptr;
  int countRowsPerProc;
  int *buffer = nullptr;
  double startTime;
  double endTime; 
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &countProc);
  startTime = MPI_Wtime();
  if (rank == 0) {
	lenRow = 10;
	matrix = new int[countRow * lenRow];
	std::mt19937 mt(1);
	for (size_t i = 0; i < countRow * lenRow; i++)
	{
	  matrix[i] = mt() % 10;
	}

	if (countRow % countProc == 0) {
	  countRowsPerProc = countRow / countProc;
	  remainder = 0;
	}
	else {
	  countRowsPerProc = countRow / countProc;
	  remainder = countRow % countProc;
	}
	result = new int[countRowsPerProc * countProc + remainder];
  }
  MPI_Bcast(&lenRow, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&countRowsPerProc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&remainder, 1, MPI_INT, 0, MPI_COMM_WORLD);
  localResult = new int[countRowsPerProc + remainder];
  buffer = new int[countRowsPerProc * lenRow + remainder * lenRow];
  MPI_Scatter(matrix, countRowsPerProc * lenRow, MPI_INT, buffer, countRowsPerProc * lenRow,
	MPI_INT, 0, MPI_COMM_WORLD);

  for (size_t i = 0; i < countRowsPerProc; i++)
  {
	localResult[i] = 0;
	for (size_t j = i * lenRow; j < lenRow + i * lenRow; j++)
	{
	  localResult[i] += buffer[j];
	}
  }
 
  MPI_Gather(localResult, countRowsPerProc, MPI_INT, result, countRowsPerProc, MPI_INT, 0,
	  MPI_COMM_WORLD);
  if (remainder != 0) {
	if (rank == 0) {
	  for (size_t i = 1; i < remainder; i++)
	  {
		MPI_Send(&matrix[countRowsPerProc * countProc * lenRow + i * lenRow], lenRow, MPI_INT, i, 1, MPI_COMM_WORLD);
	  }
	  result[countRowsPerProc * countProc] = 0;
	  for (size_t i = 0; i < lenRow; i++)
	  {
		result[countRowsPerProc * countProc] += matrix[countRowsPerProc * countProc * lenRow
		  + i];
	  }

	  for (size_t i = 1; i < remainder; i++)
	  {
		MPI_Recv(&result[countRowsPerProc * countProc + i], 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
	  }
	}
	else {
	  if (rank < remainder) {
		MPI_Recv(buffer, lenRow, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		for (size_t i = 0; i < lenRow; i++)
		{
		  localResult[countRowsPerProc * lenRow + 1] += buffer[i];
		}
		MPI_Send(&localResult[countRowsPerProc * lenRow + 1], 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	  }
	}
  }
 
  endTime = MPI_Wtime();
  if (rank == 0) {
	for (size_t i = 0; i < countRow * lenRow; i++)
	{
	  std::cout << matrix[i] << " ";
	  if (i % (lenRow - 1) == 0 && i != 0)
		std::cout << std::endl;
	}
	std::cout << std::endl;
	for (size_t i = 0; i < countRowsPerProc * countProc + remainder; i++)
	{
	  std::cout << result[i] << std::endl;
	}
	delete[] matrix;
	std::cout << "Time: " << endTime - startTime << std::endl;
  }
  MPI_Finalize();
  delete[] localResult;
  delete[] result;
}