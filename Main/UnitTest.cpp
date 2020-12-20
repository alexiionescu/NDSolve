#include <stdio.h>
#include <stdlib.h>
#include "Function.h"
using namespace FunctionalMath;
using namespace std;
#include "Plot.h"
#include "UnitTest.h"

//#define FORMULA(t)	0.3+5*exp(-t*t/2)*sin(3*t)-t+4*cos(t)+atan(t)

void UnitTest::OperatorsTest() {
	int tickn, tick = os_ticks();
	char title1[128]{}, title2[128]{};
	const char* titles[5] = { "Operators Test", title1, title2, NULL, NULL };
	
	{//real functions 
		double T = 1;
		FunctionPtr  f1(new Function(T, [T](double t)->double{
			return gauss_sin(0.1, -0.5, t, π);
		}, -T,5));
		FunctionPtr  f2(new Function(T, [T](double t)->double{
			return gauss_cos(0.1, 0.5, t, π);
		}, -T,5));

		Operator O{ [](double x, Function& f, long long s){
			return (1+pow(x,2)) * f[s];
		} };


		FunctionPtr P = O | 2-(f1+f2)/2 +(1+f1)*(1+f2);
		tickn = os_ticks();
		printf("Compute in %d ms\n", tickn - tick);
		tick = tickn;

		
		{
			Function*  funcToPlot[] = { f1.get(), f2.get(), P.get() };
			const char* legends[] = { "f1", "f2", "Expr", NULL };
			AnalyzesResult a{};
			for (int i = 0; legends[i]; i++)
			{
				printf("\n------ %s -------\n", legends[i]);
				funcToPlot[i]->Analyze(a, stdout);
			}
			Plot(sizeof(funcToPlot) / sizeof(funcToPlot[0]), funcToPlot,
				titles, "s", "", legends);
		}
	}

	
	{//complex functions
		tick = os_ticks();
		double σ = 0.2;
		double ωr = π / 4;
		double ωi = π / 2;
		double T = 2;
		ComplexFunctionPtr  f1(new ComplexFunction(T,
			[=](double t)->double{
			return gauss_cos(σ, -0.5, t, ωr);
		}, [=](double t)->double{
			return gauss_sin(σ, -0.5, t, ωi);
		}, -T, 3));
		ComplexFunctionPtr  f2(new ComplexFunction(T,
			[=](double t)->double{
			return gauss_cos(σ, 0.5, t, ωr);
		},
			[=](double t)->double{
			return gauss_sin(σ, 0.5, t, ωi);
		}, -T, 3));


		ComplexOperator<> Â{
			//real operator
			[](double x, ComplexFunction& f, long long s){
				return x*(f.real()[s]);
			},
				//imag operator
				[](double x, ComplexFunction& f, long long s){
				return x*(f.imag()[s]);
			}
		};

		Complex a(2, 0), b(-1, 1), c(0, 1);
		ComplexFunctionPtr psi = c- (a*f1 + b*f2);
		ComplexFunctionPtr P = Â | psi;  // this is Â|psi1>
		Complex expect_A_psi = psi | P;   // this is <psi1|Â|psi1>
		Complex normPsi = psi | psi;  //this should be 1 after Normalize

		tickn = os_ticks();
		printf("Compute in %d ms\n", tickn - tick);
		tick = tickn;

		titles[0] = "Operator A = x";
		snprintf(title1, sizeof(title1) - 1, 
			"$< \\psi | \\psi > =%.2lf+i(%.2lf) $", normPsi.real(), normPsi.imag());
		snprintf(title2, sizeof(title2) - 1,
			"$< \\psi |A| \\psi > =%.2lf+i(%.2lf)$", expect_A_psi.real(), expect_A_psi.imag());
		
		{
			Function*  funcToPlot[] = { psi.get()->Re(), P.get()->Re() };
			const char* legends[] = { "Re psi1", "Re A|psi1>", NULL };
			AnalyzesResult a{};
			for (int i = 0; legends[i]; i++)
			{
				printf("\n------ %s -------\n", legends[i]);
				funcToPlot[i]->Analyze(a, stdout);
			}
			Plot(sizeof(funcToPlot) / sizeof(funcToPlot[0]), funcToPlot,
				titles, "s", "", legends);
		}
		{
			Function*  funcToPlot[] = { psi.get()->Im(), P.get()->Im() };
			const char* legends[] = { "Im psi1", "Im A|psi1>", NULL };
			AnalyzesResult a{};
			for (int i = 0; legends[i]; i++)
			{
				printf("\n------ %s -------\n", legends[i]);
				funcToPlot[i]->Analyze(a, stdout);
			}
			Plot(sizeof(funcToPlot) / sizeof(funcToPlot[0]), funcToPlot,
				titles, "s", "", legends);
		}
	}
}



void UnitTest::FourierTest()
{
	double T=6;
	Function  test(T,[T](double t)->double{
		return gauss_sin(0.4, 0, t, π) + gauss_cos(0.2, -3, t, π) + gauss(0.1, 3, t);
	},-T,3);

	int tick = os_ticks();
	double MaxF=20;
	ComplexFunction * F = test.FourierTransform(MaxF,2);
	int tickn = os_ticks();
	printf("Compute Fourier in %d ms\n",tickn-tick);
	tick = tickn;

	

	const char* titles[] = {"Fourier Transform","",NULL};
	AnalyzesResult a{};

	Function * m = F->Module();
	Function *  funcToPlot[]={m,F->Re(),F->Im()};
	const char* legends[] = {"|F|","Re F","Im F",NULL};

	
	for(int i=0; legends[i];i++)
	{
		printf("\n------ %s -------\n",legends[i]);
		funcToPlot[i]->Analyze(a,stdout);
	}
	Plot(sizeof(funcToPlot)/sizeof(funcToPlot[0]),funcToPlot,titles,"Hz","",legends);
	delete m;

	Function * invF = F->InverseFourierTransform(T,-T,3);
	tickn = os_ticks();
	printf("Compute Inv Fourier in %d ms\n",tickn-tick);
	tick = tickn;
	Function *  funcToPlotInv[]={&test,invF};
	const char* legendsInv[] = {"f","inv F",NULL};

	for(int i=0; legendsInv[i];i++)
	{
		printf("\n------ %s -------\n",legendsInv[i]);
		funcToPlotInv[i]->Analyze(a,stdout);
	}
	Plot(sizeof(funcToPlotInv)/sizeof(funcToPlotInv[0]),funcToPlotInv,titles,"s","",legendsInv);
	delete invF;

	delete F;
}

void UnitTest::MultiArgTest(){
	int tickn, tick = os_ticks();
	double XMin[] {0, 0, 0};
	double XMax[] {1, 1, 1};
	int p[] {3, 3, 3};
	Function  f(3, XMax, XMin,
		[&XMin](double* X)->double {
		return  pow(X[0], 2) + pow(X[1], 2) + pow(X[2], 2);
	}, p);
	FunctionPtr dx0_f = FunctionPtr(f.CalculatePartialDerivate(0));
	FunctionPtr dx1_f = FunctionPtr(f.CalculatePartialDerivate(1));
	FunctionPtr dx2_f = FunctionPtr(f.CalculatePartialDerivate(2));
	FunctionPtr d2x1_f = FunctionPtr(f.CalculatePartialDerivate(1, 2));
	double integral = f.Integrate(XMin, XMax);
	
	tickn = os_ticks();
	printf("Compute in %d ms\n", tickn - tick);
	tick = tickn;

	AnalyzesResult a{};
	const int S = 3;
	double X[][S]{{0, 0, 0}, { 0.5, 0.2, 0.7 }, { 0.99, 0.9, 0.3 }};
	double U[]{0, 1};
	double V[]{1, 0};

	
	for (int s = 0; s < S; s++) {
		printf("\nX = {%lg, %lg, %lg}\n", X[s][0], X[s][1], X[s][2]);
		printf("f[X] = %lg\n",  f[X[s]]);
		printf("pd(x0) = %lg\n", (*dx0_f)[X[s]]);
		printf("pd(x1) = %lg\n", (*dx1_f)[X[s]]);
		printf("pd(x2) = %lg\n", (*dx2_f)[X[s]]);
		printf("pd2(x1) = %lg\n", (*d2x1_f)[X[s]]);
	}

	printf("Integral f = %lg\n", integral);


	Function& fU = f(U);
	Function& fV = f(V);
	
	const char* titles[] = { "f(X)", "", NULL };
	Function *  funcToPlot[] = {&fU , &fV };

	tickn = os_ticks();
	printf("Compute retrieve values in %d ms\n", tickn - tick);
	tick = tickn;

	char l1[100],l2[100];
	snprintf(l1, sizeof(l1) - 1, "f({%lg, %lg})", U[0], U[1]);
	snprintf(l2, sizeof(l2) - 1, "f({%lg, %lg})", V[0], V[1]);
	const char* legends[] = { l1, l2, NULL };
	for (int i = 0; legends[i]; i++)
	{
		printf("\n------ %s -------\n", legends[i]);
		funcToPlot[i]->Analyze(a, stdout);
	}
	Plot(sizeof(funcToPlot) / sizeof(funcToPlot[0]), funcToPlot, titles, "x0", "", legends);
}

/* solve 

Heat equation
	f:[0,1]x[0,1]->R

	∂f/∂x1-∂^2f/∂x0^2==0

	f(0,x1) = 1+exp(-x1)
	f(x0,0) = 2+exp(x0) 
	∂f/∂x0(x0,1) = 2*exp(1-x0)

*/

double _PDEEqn(double* X, Function::PDEData& D, void*) {
	
	return D(1) - D(0, 0);
}

void UnitTest::SolvePDE() {

	double XMin[] {0, 0, 0, 0};
	double XMax[] {0.5, 0.5, 1, 1};
	int p[] {3, 3, 3, 3};
	int const n = 2;

	Function  f(n, XMax, XMin,
		[](double* X)->double {
		if (X[0] == 0) return exp(X[1]); // IC
			
		return HUGE_VAL; //default, uninitialized
	}, p);

	int maxSolutions = 1;
	Function::DBG_PDE_SOLVE = true;
	unsigned int ticks = os_ticks();
	Function** fs = f.SolvePDENL(maxSolutions, _PDEEqn, NULL, 1, 2000, 400, 6,0);
	printf("\n --- Found %d solutions in %d ms -- \n\n", maxSolutions, os_ticks() - ticks);

	

	double v[][n - 1]{
					  {XMin[0] + 1. / Function::__Pow10(p[0])}, 
					  { (XMin[0] + XMax[0]) / 2 }, 
					  { XMax[0] - 2. / Function::__Pow10(p[0]) } 
	};
	double u[] { XMin[n - 1] + 1. / Function::__Pow10(p[n - 1]), 
				 (XMin[n - 1] + XMax[n - 1]) / 2 , 
				 XMax[n - 1] - 2. / Function::__Pow10(p[n - 1]) 
	};
	const int vsize = sizeof(v) / sizeof(v[0]);
	

	char* legends[vsize + 1]; 
	legends[vsize] = NULL;
	char* titles[] = { new char[100], "", NULL };


	for (int j = 0; j < vsize; j++) {
		legends[j] = new char[100];
		snprintf(legends[j], 99, "sol({%lg,...})", v[j][0]);
	}

	for (int k = 0; k < maxSolutions; k++) {
		
		FunctionPtr d1(fs[k]->CalculatePartialDerivate(1));
		FunctionPtr d00(fs[k]->CalculatePartialDerivate(0, 2));

		//print some values
		printf("\n");
		int rvsize = vsize;
		for (int j = 0; j < rvsize; j++) {
			for (int l = 0; l < sizeof(u) / sizeof(u[0]); l++) {
				double val = (*fs[k])(v[j])[u[l]];
				if (val == HUGE_VAL){
					rvsize = j;
					break;
				}
				printf("sol%d({ %lg,...,%lg }) = %lg\n", k,v[j][0], u[l], val);
				printf("d1 sol%d({ %lg,...,%lg }) = %lg\n", k,v[j][0], u[l], (*d1)(v[j])[u[l]]);
				printf("d00 sol%d({ %lg,...,%lg }) = %lg\n\n", k,v[j][0], u[l], (*d00)(v[j])[u[l]]);
			}
		}

		snprintf(titles[0], 99, "sol %d", k);

		Function *  funcToPlot[2*vsize];
		for (int j = 0; j < rvsize; j++) {
			funcToPlot[j] = &(*fs[k])(v[j]);
		}

		//for (int i = 0; legends[i]; i++)
		//{
		//	printf("\n------ %s -------\n", legends[i]);
		//	AnalyzesResult a{};
		//	funcToPlot[i]->Analyze(a, stdout);
		//}
			
		Plot(rvsize, funcToPlot, (const char**)titles, "x1", "", (const char**)legends);
		delete fs[k];
	}

	for (int j = 0; j < vsize; j++) {
		delete[] legends[j];
	}

	delete titles[0];

	delete[] fs;
}