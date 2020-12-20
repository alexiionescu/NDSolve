#pragma once

#include "Function.h"
class QuantumMechanic
{
	
public:
	static int const MAX_EIGENVALUES = 24;
	
	FunctionalMath::ComplexFunctionPtr _Ψ[MAX_EIGENVALUES * 2];
	int _precision = 4;
	//the n-th eigenvalue is requested to calculate and plot, use -1 to calculate all
	void BoundStates(const int n=-1, bool SAVE_ALL_GRAPHS = false);
	void OperatorsTest(FunctionalMath::ComplexFunctionPtr&);

	void BoundStateAndOperators();
};

