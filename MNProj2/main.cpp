#include <iostream>
#include "Matrix.h"
#include <cmath>

using namespace std;

//index: 188776
const int c = 7;
const int d = 6;
const int e = 7;
const int f = 8;
const int N = c * d;

void zadanieB(Matrix& A, Matrix& b)
{
	cout << "GaussSeidel:\n";
	EquationsSystemSolution solutionGauss = Matrix::GaussSeidel(A, b, 10e-9, 1000);
	cout << solutionGauss << endl;

	cout << "Jacobi:\n";
	EquationsSystemSolution solutionJacobi = Matrix::Jacobi(A, b, 10e-9, 1000);
	cout << solutionGauss << endl;

	// Wykresy
}

void zadanieC()
{
	Matrix A = Matrix::zeros(N, N);
	for (int row = 0; row < N; row++)
	{
		if (row - 1 >= 0)
			A.at(row, row - 1) = -1;
		if (row - 2 >= 0)
			A.at(row, row - 2) = -1;

		A.at(row, row) =3;

		if (row + 1 < N)
			A.at(row, row + 1) = -1;
		if (row + 2 < N)
			A.at(row, row + 2) = -1;
	}

	Matrix b = Matrix::zeros(N, 1);
	for (int row = 0; row < N; row++)
	{
		b.at(row, 0) = sin(row * f + 1);
	}

	EquationsSystemSolution solutionGauss = Matrix::GaussSeidel(A, b, 10e-9, 1000);
	EquationsSystemSolution solutionJacobi = Matrix::Jacobi(A, b, 10e-9, 1000);
	// Wykresy
}

void zadanieD()
{
	Matrix A = Matrix::zeros(N, N);
	for (int row = 0; row < N; row++)
	{
		if (row - 1 >= 0)
			A.at(row, row - 1) = -1;
		if (row - 2 >= 0)
			A.at(row, row - 2) = -1;

		A.at(row, row) = 5+e;

		if (row + 1 < N)
			A.at(row, row + 1) = -1;
		if (row + 2 < N)
			A.at(row, row + 2) = -1;
	}

	Matrix b = Matrix::zeros(N, 1);
	for (int row = 0; row < N; row++)
	{
		b.at(row, 0) = sin(row * f + 1);
	}

	EquationsSystemSolution solutionLU = Matrix::LU(A, b);

	// Wykresy
}

void zadanieE()
{
	const int N_length = 3;
	int N[] = { 100, 500, 1000};
	EquationsSystemSolution solutionsGauss[N_length];
	EquationsSystemSolution solutionsJacobi[N_length];


	for (int i = 0; i < N_length; i++)
	{
		Matrix A = Matrix::zeros(N[i], N[i]);
		for (int row = 0; row < N[i]; row++)
		{
			if (row - 1 >= 0)
				A.at(row, row - 1) = -1;
			if (row - 2 >= 0)
				A.at(row, row - 2) = -1;
			A.at(row, row) = 5 + e;
			if (row + 1 < N[i])
				A.at(row, row + 1) = -1;
			if (row + 2 < N[i])
				A.at(row, row + 2) = -1;
		}
		Matrix b = Matrix::zeros(N[i], 1);
		for (int row = 0; row < N[i]; row++)
			b.at(row, 0) = sin(row * f + 1);

		solutionsGauss[i] = Matrix::GaussSeidel(A, b, 10e-9, 1000);
		solutionsJacobi[i] = Matrix::Jacobi(A, b, 10e-9, 1000);
	}
	cout << "Gauss: ";
	for(int i = 0; i < N_length; i++)
		cout << N[i] << ":\n\t" << solutionsGauss[i] << endl;

}

int main()
{

	cout.precision(2);
	
	// Zadanie A
	Matrix A = Matrix::zeros(N, N);
	for (int row = 0; row < N; row++)
	{
		if (row - 1 >= 0)
			A.at(row, row - 1) = -1;

		A.at(row, row) = 5 + e;

		if (row + 1 < N)
			A.at(row, row + 1) = -1;
	}

	Matrix b = Matrix::zeros(N, 1);
	for (int row = 0; row < N; row++)
	{
		b.at(row, 0) = sin(row * f + 1);
	}

	//zadanieB(A, b);
	//zadanieC();
	//zadanieD();
	zadanieE();
	
}