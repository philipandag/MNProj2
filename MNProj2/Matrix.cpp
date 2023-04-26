#include "Matrix.h"

template <typename T>
T& min(T& a, T& b){
	return a < b ? a : b;
}
template <typename T>
T& max(T& a, T& b) {
	return a > b ? a : b;
}
template <typename T>
void swap(T& a, T& b)
{
	T temp = a;
	a = b;
	b = temp;
}

Matrix::Matrix():
	Matrix(0, 0)
{}

Matrix::Matrix(const Matrix& other) :
	Matrix(other.rows, other.columns, other.data)
{
	for (int i = 0; i < rows * columns; i++)
	{
		data[i] = other.data[i];
	}
}

Matrix::Matrix(Matrix&& other) :
	rows(other.rows),
	columns(other.columns),
	data(std::move(other.data))
{}

Matrix::Matrix(int rows, int columns):
	rows(rows),
	columns(columns),
	data(new double[rows*columns])
{}

Matrix::Matrix(int N) :
	Matrix(N, N)
{}

Matrix::Matrix(int rows, int columns, double* data):
	rows(rows),
	columns(columns),
	data(data)
{}

Matrix& Matrix::operator=(const Matrix& other) {
	rows = other.rows;
	columns = other.columns;
	data = new double[rows * columns];
	for (int i = 0; i < rows * columns; i++)
		data[i] = other.data[i];
	
	return *this;
}
/*
Matrix& Matrix::operator=(Matrix other)
{
	std::swap(rows, other.rows);
	std::swap(columns, other.columns);
	std::swap(data, other.data);
	return *this;
}
*/
Matrix& Matrix::operator=(Matrix&& other) {
	rows = other.rows;
	columns = other.columns;
	data = std::move(other.data);
	return *this;
}


double& Matrix::at(int row, int column)
{
	return data[row * columns + column];
}

const double& Matrix::at(int row, int column) const
{
	return data[row * columns + column];
}

Matrix Matrix::copy() const
{
	Matrix m(rows, columns);
	for (int i = 0; i < rows * columns; i++)
	{
		m.data[i] = data[i];
	}
	return m;
}

Matrix& Matrix::operator+=(Matrix& other)
{
	if(rows!=other.rows || columns != other.columns)
		throw "Matrix adding, dimensions dont agree";
	for (int i = 0; i < rows * columns; i++)
	{
		data[i] += other.data[i];
	}
	return *this;
}

Matrix Matrix::add(const Matrix& one, const Matrix& other)
{
	if (one.rows != other.rows || one.columns != other.columns)
		throw "Matrix adding, dimensions dont agree";

	Matrix result(one.rows, one.columns);
	for (int i = 0; i < one.rows * one.columns; i++)
	{
		result.data[i] = one.data[i] + other.data[i];
	}
	return result;
}

Matrix Matrix::subtract(const Matrix& one, const Matrix& other)
{
	int minRows = min(one.rows, other.rows);
	int minColumns = min(one.columns, other.columns);

	Matrix result(minRows, minColumns);

	for (int row = 0; row < minRows; row++)
		for (int column = 0; column < minColumns; column++)
			result.at(row, column) = one.at(row, column) - other.at(row, column);
	
	return result;
}

Matrix Matrix::subtract(const Matrix& one, double value)
{
	Matrix result(one.rows, one.columns);
	for (int row = 0; row < one.rows; row++)
		for (int column = 0; column < one.columns; column++)
			result.at(row, column) = one.at(row, column) - value;

	return result;
}

Matrix Matrix::dotProduct(const Matrix& one, const  Matrix& other)
{
	if (one.columns != other.rows)
		throw "Matrix dot, dimensions dont agree";
	Matrix result(one.rows, other.columns);
	for (int row = 0; row < one.rows; row++)
	{
		for (int column = 0; column < other.columns; column++)
		{
			double sum = 0;
			for (int i = 0; i < one.columns; i++)
			{
				sum += one.at(row, i) * other.at(i, column);
			}
			result.at(row, column) = sum;
		}
	}
	return result;
}

Matrix& Matrix::T()
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			swap(at(i, j), at(j, i));
		}
	}
	return *this;
}
Matrix Matrix::transpose(const Matrix& m)
{
	Matrix result(m.columns, m.rows);
	for (int i = 0; i < m.rows; i++)
	{
		for (int j = 0; j < m.columns; j++)
		{
			result.at(i, j) = m.at(j, i);
		}
	}
	return result;
}

Matrix Matrix::eye(int n)
{
	
	Matrix m(n);
	for (int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			m.at(i, j) = i == j ? 1 : 0;
	return m;
}

Matrix Matrix::zeros(int rows, int columns)
{
	Matrix m(rows, columns);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < columns; j++)
			m.at(i, j) = 0;
	return m;
}

Matrix Matrix::ones(int rows, int columns)
{
	Matrix m(rows, columns);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < columns; j++)
			m.at(i, j) = 1;
	return m;
}

Matrix Matrix::forwardSubstitution(const Matrix& A, const Matrix& b)
{
	if (A.rows != A.columns || A.rows != b.rows)
		throw "Matrix forwardSubstitution, dimensions dont agree";

	Matrix x(b.rows, 1);
	for (int i = 0; i < b.rows; i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
		{
			sum += A.at(i, j) * x.at(j, 0);
		}
		x.at(i, 0) = (b.at(i, 0) - sum) / A.at(i, i);
	}
	return x;
}

Matrix Matrix::backwardSubstitution(const Matrix& A, const Matrix& b)
{
	if (A.rows != A.columns || A.rows != b.rows)
		throw "Matrix backwardSubstitution, dimensions dont agree";
	Matrix x(b.rows, 1);
	for (int i = b.rows - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < b.rows; j++)
		{
			sum += A.at(i, j) * x.at(j, 0);
		}
		x.at(i, 0) = (b.at(i, 0) - sum) / A.at(i, i);
	}
	return x;
}

Matrix Matrix::upperTriangular(const Matrix& A, int offsetFromDiagonal)
{
	if (A.rows != A.columns)
		throw "Matrix upperTriangular, dimensions dont agree";

	Matrix result(A.rows, A.columns);

	for (int row = 0; row < A.rows-offsetFromDiagonal; row++)
	{
		for (int col = row + offsetFromDiagonal; col < A.columns; col++)
		{
			result.at(row, col) = A.at(row,col);
		}
	}
	return result;
}

Matrix Matrix::lowerTriangular(const Matrix& A, int offsetFromDiagonal)
{
	if (A.rows != A.columns)
		throw "Matrix lowerTriangular, dimensions dont agree";
	Matrix result(Matrix::zeros(A.rows, A.columns));

	for (int row = offsetFromDiagonal; row < A.rows; row++)
	{
		for (int col = 0; col < row + 1 - offsetFromDiagonal; col++)
		{
			result.at(row, col) = A.at(row, col);
		}
	}
	return result;
}

Matrix Matrix::diagonal(const Matrix& A)
{
	if (A.rows != A.columns)
		throw "Matrix diagonal, dimensions dont agree";

	Matrix result(Matrix::zeros(A.rows, A.columns));
	for (int row = 0; row < A.rows; row++)
	{
		result.at(row, row) = A.at(row, row);
	}
	return result;
}

double Matrix::residualError(const Matrix& A, const Matrix& b, const Matrix& x)
{
	Matrix Ax(Matrix::dotProduct(A, x));
	Matrix r(Matrix::subtract(b, Ax));

	double sum = 0;
	for (int i = 0; i < r.rows; i++)
	{
		sum += r.at(i, 0) * r.at(i, 0);
	}

	return sqrt(sum);
}

EquationsSystemSolution Matrix::GaussSeidel(const Matrix& A, const Matrix& b, double tolerance, int maxIterations)
{
	clock_t t = clock();

	Matrix DL(Matrix::lowerTriangular(A));
	Matrix U(Matrix::upperTriangular(A,1));
	Matrix x(Matrix::ones(b.rows, 1));

	double* residualErrors = new double[maxIterations];
	int iterations = 0;

	for (int i = 0; i < maxIterations; i++)
	{
		iterations += 1;

		Matrix Ux(Matrix::dotProduct(U, x));
		Matrix bmUx(Matrix::subtract(b, Ux));

		x = Matrix::forwardSubstitution(DL, bmUx);


		double error = Matrix::residualError(A, b, x);
		residualErrors[i] = error;
		if (Matrix::residualError(A, b, x) < tolerance)
			break;
	}

	return EquationsSystemSolution(iterations, residualErrors,
		(clock() - t) / (double)CLOCKS_PER_SEC, &x);
}

EquationsSystemSolution Matrix::Jacobi(const Matrix& A, const Matrix& b, double tolerance, int maxIterations)
{
	clock_t t = clock();

	Matrix D(Matrix::diagonal(A));
	Matrix LU(Matrix::add(Matrix::lowerTriangular(A, 1), Matrix::upperTriangular(A, 1)));

	Matrix x = Matrix::ones(b.rows, 1);

	double* residualErrors = new double[maxIterations];
	int iterations = 0;

	for (int i = 0; i < maxIterations; i++)
	{
		iterations += 1;
		Matrix LUx(Matrix::dotProduct(LU, x));
		Matrix bmLUx(Matrix::subtract(b, LUx));
		Matrix newX(Matrix::forwardSubstitution(D, bmLUx));
		
		x = newX;

		double error = Matrix::residualError(A, b, x);
		residualErrors[i] = error;
		if (Matrix::residualError(A, b, x) < tolerance)
			break;
	}

	return EquationsSystemSolution(iterations, residualErrors,
		(clock() - t) / (double)CLOCKS_PER_SEC, &x);
}

void Matrix::LUDecomposition(const Matrix& A, Matrix* L, Matrix* U)
{
	if (A.rows != A.columns)
		throw "Matrix LUDecomposition, dimensions dont agree";
	*L = Matrix::eye(A.rows);
	*U = A.copy();


	int m = A.rows;
	for (int k = 0; k < m - 1; k++)
	{
		for (int j = k + 1; j < m; j++)
		{
			L->at(j, k) = U->at(j, k) / U->at(k, k);
			for (int i = k; i < m; i++)
			{
				U->at(j, i) = U->at(j, i)- L->at(j, k) * U->at(k, i);
			}
		}
	}
	return;
}

EquationsSystemSolution Matrix::LU(const Matrix& A, const Matrix& b)
{
	clock_t t0 = clock();
	Matrix L;
	Matrix U;
	Matrix::LUDecomposition(A, &L, &U);
	Matrix y(Matrix::forwardSubstitution(L, b));
	Matrix x(Matrix::backwardSubstitution(U, y));

	return EquationsSystemSolution(-1, nullptr, 
		(clock() - t0) / (double)(CLOCKS_PER_SEC), &x);
}

Matrix::~Matrix()
{
	if(data != nullptr)
		delete[] data;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m)
{
	for (int i = 0; i < m.rows; i++)
	{
		for (int j = 0; j < m.columns; j++)
		{
			os << std::scientific << m.at(i, j) << " ";
		}
		os << std::endl;
	}
	return os;
}
