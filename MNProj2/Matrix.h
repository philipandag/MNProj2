#pragma once
#include "IterationMethodResult.h"
#include <iostream>
#include <time.h>

struct Matrix
{
	int rows, columns;
	double* data;

	Matrix();
	Matrix(const Matrix& other);
	Matrix(Matrix&& other);
	Matrix(int rows, int columns);
	Matrix(int N);
	Matrix(int rows, int columns, double* data);

	Matrix& operator=(const Matrix& other);
	Matrix& operator=(Matrix&& other);
	//Matrix& operator=(Matrix other);


	static Matrix eye(int n);
	static Matrix zeros(int rows, int columns);
	static Matrix ones(int rows, int columns);
	double& at(int row, int column);
	const double& at(int row, int column) const;

	//creates new matrix
	Matrix copy() const;

	//modifies this matrix
	Matrix& operator+=(Matrix& other);
	//creates new matrix
	static Matrix add(const Matrix& one, const Matrix& other);
	static Matrix subtract(const Matrix& one, const Matrix& other);
	static Matrix subtract(const Matrix& one, const double value);
	
	//modifies this matrix
	Matrix& T();

	//creates new matrix
	static Matrix transpose(const Matrix& m);

	// creates new matrix
	static Matrix dotProduct(const Matrix& one, const Matrix& other);

	// creates new matrix
	static Matrix forwardSubstitution(const Matrix& A, const Matrix& b);
	static Matrix backwardSubstitution(const Matrix& A, const Matrix& b);

	// creates new matrix
	static Matrix upperTriangular(const Matrix& A, int offsetFromDiagonal=0);
	static Matrix lowerTriangular(const Matrix& a, int offsetFromDiagonal=0);
	static Matrix diagonal(const Matrix& A  );

	static double residualError(const Matrix& A, const Matrix& b, const Matrix& x);

	static EquationsSystemSolution GaussSeidel(const Matrix& A, const Matrix& b, double tolerance = 0.0001, int maxIterations = 1000);
	static EquationsSystemSolution Jacobi(const Matrix& A, const Matrix& b, double tolerance = 0.0001, int maxIterations = 1000);

	static void LUDecomposition(const Matrix& A, Matrix* L, Matrix* U);
	static EquationsSystemSolution LU(const Matrix& A, const Matrix& b);
	~Matrix();

	friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
};

