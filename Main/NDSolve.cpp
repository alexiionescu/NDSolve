// NDSolve.cpp : Defines the entry point for the console application.
//
#include <stdio.h>
#include "n-body-problem.h"
#include "QuantumMechanic.h"
#include "UnitTest.h"
#include "Function.h"


int main(int argc, char* argv[])
{
	printf("NDSolve IN with Function::maxN=%I64u\n", FunctionalMath::Function::maxN);
	try{
		UnitTest ut;
		ut.OperatorsTest();
		//ut.MultiArgTest();
		//ut.FourierTest();
		//ut.SolvePDE(); //TODO

		//NBodyProblem nbp;
		//nbp.NewObjectSolarSystem();

		//QuantumMechanic qm;
		//qm.BoundStateAndOperators();	
	}
	catch (const char* err){
		printf("\nEXCEPTION %s\n", err);
	}
	catch (...){
		printf("\nUNKNOWN EXCEPTION\n");
	}
	printf("\nNDSolve OUT with Function::maxN=%I64u\n", FunctionalMath::Function::maxN);
	return 0;
}

