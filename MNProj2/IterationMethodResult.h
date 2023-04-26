#pragma once
#include <time.h>
#include <iostream>

struct Matrix;

struct EquationsSystemSolution
{
	int iterations;
	double* errorArray;
	double time;
	Matrix* solution;
	EquationsSystemSolution();
	EquationsSystemSolution(int iterations, double* errorArray, double time, Matrix* solution);
	~EquationsSystemSolution();
	friend std::ostream& operator<<(std::ostream& os, const EquationsSystemSolution& m);
};