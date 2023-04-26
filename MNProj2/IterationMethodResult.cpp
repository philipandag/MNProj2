#include "IterationMethodResult.h"
#include "Matrix.h"

EquationsSystemSolution::EquationsSystemSolution():
	iterations(-1),
	errorArray(nullptr),
	time(-1),
	solution(nullptr)
{}

EquationsSystemSolution::EquationsSystemSolution(int iterations, double* errorArray, double time, Matrix* solution):
	iterations(iterations),
	errorArray(errorArray),
	time(time),
	solution(new Matrix(solution->copy()))
{
	solution->data = nullptr;
}

EquationsSystemSolution::~EquationsSystemSolution()
{
	if(errorArray != nullptr)
		delete [] errorArray;
	if(solution != nullptr)
		delete solution;
}

std::ostream& operator<<(std::ostream& os, const EquationsSystemSolution& m)
{
	os << "Iterations: " << m.iterations << std::endl;
	os << "Time: " << m.time << std::endl;
	os << "Error: " << m.errorArray[m.iterations-1] << std::endl;
	return os;
}