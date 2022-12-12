#include <iostream>
#include "mpi.h"
#include <vector>
//int main(int argc, char* argv[]) {
//
//  double* localB; // Те составляющие, которые необходимы для вычисления(приходят 1 раз)
//  std::vector<double> equations; // Те кэфы при составляющей х, которую нужно посчитать(строки матрицы А)(1 раз)
//  //std::vector<double> predX;//Вычисленные на прошлых итерациях составляющие х
//  //std::vector<double> localX;
//  //double* equations;
//  double* predX = nullptr;
//  std::vector<double> diag;
//  double* localX;
//  const int size = 3; // передается откуда-то из вне, как и матрица. - это размер матрицы A.
//  const int lenRow = size + 1;
//  int rank;
//  diag.resize(size);
//  int countProc;
//  int remainder; // Для понимания того, через scatter или scatter_v будем действовать
//  int dataPerProc;
//  //std::vector<int> sendCounts;
//  //std::vector<int> displs;
//  int* sendCounts = nullptr;
//  int* displs = nullptr;
//  int mpi_recv;
//  double accuracy = 0.000001;
//  double norm;
//  //double* matrix = new double[size * lenRow]; //Где-то внутри вызывается метод GetRandomMatrix()
//  double matrix[] = { 4.0, -2.0, 1.0, 0,
//					2.0, 5.0, -1.0, 3.0,
//					4.0, -1.0, 7.0, 5.0 };
//
//  //Блок инициализации
//  MPI_Init(&argc, &argv);
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  MPI_Comm_size(MPI_COMM_WORLD, &countProc);
//  dataPerProc = size / countProc;
//  remainder = size % countProc;
//
//  if (remainder != 0) {
//	if (rank < remainder) {
//	  mpi_recv = dataPerProc + 1;
//	}
//	else {
//	  mpi_recv = dataPerProc;
//	}
//  }
//  else {
//	mpi_recv = dataPerProc;
//  }
//
//  //выделения места для начального приближени
//  predX = new double[size];
//
//
//
//  // Цели: узнать, что используем s or s_v.
//  // Передать каждому процессу строки матрицы(equations)
//  // Передать каждому процессу составляющие вектора b(сделать за один шаг)
//  // Задать первое приближение(заполнить массив predX) и передать другим процессам(как?)
//  // (скорее всего Bcast)
//
//
//  //нужно где-то отдельно обрабатывать ситуацию, когда число процессов больше числа строк
//
//  // начальная работа с передачей
// // if (rank == 0) {
//	//
//
//	//if (remainder == 0) {
//	//  //передача строк матрицы
//	//  equations = new double[lenRow * dataPerProc];
//	//  MPI_Scatter(matrix, dataPerProc * lenRow, MPI_INT, equations, dataPerProc * lenRow,
//	//	MPI_INT, 0, MPI_COMM_WORLD);
//	//  //работа с b
//	//  localB = new double[dataPerProc];
//	//  for (size_t i = 0; i < dataPerProc; i++)
//	//	localB[i] = equations[size];
//
//
//
//	//  //К этому моменту передано: начальне приближение, строки матрицы A|b
//	//}
//	//else {
//	//  sendCounts = new int[countProc];
//	//  displs = new int[countProc];
//	//  sendCounts[0] = lenRow * (dataPerProc + 1);
//	//  displs[0] = 0;
//	//  for (size_t i = 1; i < countProc; i++)
//	//  {
//	//	if (rank < remainder) {
//	//	  sendCounts[i] = (dataPerProc + 1) * lenRow;
//	//	  displs[i] = rank * (dataPerProc + 1) * lenRow;
//	//	}
//	//	else {
//	//	  sendCounts[i] = dataPerProc * lenRow;
//	//	  displs[i] = rank * dataPerProc + remainder;
//	//	}
//	//  }
//	//}
//	//MPI_Scatterv(matrix, sendCounts, displs, MPI_INT, equations, sendCounts[rank], MPI_INT, 0,
//	//  MPI_COMM_WORLD);
//	//// Получение вектора b
//	//for (size_t i = 0; i < dataPerProc; i++)
//	//  localB[i] = equations[size];
// // }
// // else {
//	//MPI_Bcast(&predX, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	//if (remainder == 0) {
//	//  equations = new double[lenRow * dataPerProc];
//	//  MPI_Scatter(matrix, dataPerProc * lenRow, MPI_INT, equations, dataPerProc * lenRow,
//	//	MPI_INT, 0, MPI_COMM_WORLD);
//	//  //получение вектора b
//	//  for (size_t i = 0; i < dataPerProc; i++)
//	//	localB[i] = equations[size];
//
//	//}
//	//else {
//	//  
//	//  sendCounts[0] = lenRow * (dataPerProc + 1);
//	//  displs[0] = 0;
//	//  for (size_t i = 1; i < countProc; i++)
//	//  {
//	//	if (rank < remainder) {
//	//	  sendCounts[i] = (dataPerProc + 1) * lenRow;
//	//	  displs[i] = rank * (dataPerProc + 1) * lenRow;
//	//	}
//	//	else {
//	//	  sendCounts[i] = dataPerProc * lenRow;
//	//	  displs[i] = rank * dataPerProc + remainder;
//	//	}
//	//  }
//	//  MPI_Scatterv(matrix, sendCounts, displs, MPI_INT, equations, sendCounts[rank], MPI_INT, 0,
//	//	MPI_COMM_WORLD);
//	//}
//	//// Получение вектора b
//	//
// // }
//  if (rank == 0) {
//	//задание начального приближения
//	for (size_t i = 0; i < size; i++) {
//	  predX[i] = 0;//matrix[i * lenRow + size];
//	}	
//	// передача начального приближения всем процессам
//	MPI_Bcast(predX, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	sendCounts = new int[countProc];
//	displs = new int[countProc];
//	int r = remainder;
//	for (int iRank = 0; iRank < countProc; iRank++) {
//	  if (r != 0) {
//		sendCounts[iRank] = (dataPerProc + 1) * lenRow;
//		displs[iRank] = iRank * ((dataPerProc + 1) * lenRow);
//		r--;
//	  }
//	  else {
//		sendCounts[iRank] = dataPerProc * lenRow;
//		displs[iRank] = (iRank + remainder) * dataPerProc * lenRow;
//	  }
//	}
//	
//  }
//  else {
//	MPI_Bcast(predX, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//  }
// 
//  localX = new double[mpi_recv];
//  // Передача строк матрицы A|b
//  //equations = new double[mpi_recv * lenRow];
//  equations.resize(mpi_recv * lenRow);
//  MPI_Scatterv(matrix, sendCounts, displs, MPI_DOUBLE, &equations[0], mpi_recv * lenRow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//  /*if (rank == 0) {
//	for (size_t i = 0; i < mpi_recv * lenRow; i++)
//	{
//	  std::cout << "My rank:" << 0 << std::endl;
//	  std::cout << equations[i] << " ";
//	}
//  }*/
//  if (rank == 0) {
//	std::cout << "Before:" << std::endl;
//	std::cout << "My rank:" << 0 << std::endl;
//	for (size_t i = 0; i < mpi_recv * lenRow; i++)
//	{
//	  std::cout << equations[i] << " ";
//	}
//  }
//
//
//
//  //Получение вектора b
//  /*for (size_t i = 0; i < mpi_recv; i++)
//	localB[i] = equations[i * lenRow + size];*/
//   //преобразование коэффициентов для того, чтобы потом только изменять х
//
//
//  for (size_t row = 0; row < mpi_recv; row++) {
//	for (size_t column = row * lenRow; column < row * lenRow + lenRow - 1; column++)
//	{
//	  if (rank < remainder) {
//		if (column == rank * (dataPerProc + 1) + row * (lenRow + 1))
//		  continue;
//		equations[column] /= equations[rank * (dataPerProc + 1) + row * (lenRow + 1)];
//		equations[column] *= -1;
//	  }
//	  else {
//		if (column == rank * (dataPerProc + remainder) + row * (lenRow + 1))
//		  continue;
//		equations[column] /= equations[rank * (dataPerProc + remainder) + row * (lenRow + 1)];
//		equations[column] *= -1;
//	  }
//	}
//	if (rank < remainder) {
//	  equations[row * lenRow + size] /= equations[rank * (dataPerProc + 1) + row * (lenRow + 1)];
//	  //занулили диагональный элемент
//	  equations[rank * (dataPerProc + 1) + row * (lenRow + 1)] = 0;
//	}
//	else {
//	  equations[row * lenRow + size] /= equations[rank * (dataPerProc + remainder) + row * (lenRow + 1)];
//	  //занулили диагональный элемент
//	  equations[rank * (dataPerProc + remainder) + row * (lenRow + 1)] = 0;
//	}
//
//  }
//
//
// 
//
//  if (rank == 0) {
//	std::cout << "After:" << std::endl;
//	for (size_t row = 0; row < mpi_recv * lenRow; row++) {
//	  std::cout << equations[row] << " ";
//	}
//  }
//
//  int k = 8;
//  // Сам метод
//
//  int* recvCounts = new int[countProc];
//  int* recvDispls = new int[countProc];
//  for (size_t i = 0; i < countProc; i++)
//  {
//	if (rank < remainder) {
//	  recvCounts[i] = dataPerProc + 1;
//	  recvDispls[i] = i * (dataPerProc + 1);
//	} 
//	else {
//	  recvCounts[i] = dataPerProc;
//	  recvDispls[i] = i * (dataPerProc + remainder);
//	}
//  }
//  
//
//  do {
//	for (size_t x_i = 0; x_i < mpi_recv; x_i++)
//	{
//	  localX[x_i] = 0;
//	  for (size_t column = x_i * lenRow; column < x_i * lenRow + lenRow; column++)
//	  {
//
//		localX[x_i] += equations[column] * predX[column % lenRow];
//	  }
//	  localX[x_i] += equations[x_i * lenRow + size];
//	}
//	if (rank == 0) {
//	  std::cout << "Iteration :" << k << std::endl;
//	  for (size_t i = 0; i < size; i++)
//	  {
//		std::cout<< predX[i] << std::endl;
//	  }
//	}
//	MPI_Allgatherv(localX, mpi_recv, MPI_DOUBLE, predX, recvCounts, recvDispls, MPI_DOUBLE, MPI_COMM_WORLD);
//
//	k--;
//  } while (k > 0);
//
//  if (rank == 0) {
//	for (size_t i = 0; i < size; i++)
//	{
//	  std::cout << predX[i] << std::endl;
//	}
//  }
//
//  MPI_Finalize();
//
//}
