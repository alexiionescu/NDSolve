#include <stdio.h>
#include <stdlib.h>
#include "UnitTest.h"
#include "Plot.h"
#include "QuantumMechanic.h"
using namespace FunctionalMath;

//environment constants
static double m = 1; //mass of particle in multiple of electron mass
static double X0 = 0;
static double X = 5; //compute length from X0 to X in nm

void QuantumMechanic::BoundStateAndOperators(){
	ComplexFunctionPtr Ψ;

	Complex C[MAX_EIGENVALUES * 2] { 
		{1, 0}, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
	};

	BoundStates();
	if (_Ψ[0].get() != NULL) {
		Ψ = C[0] * _Ψ[0];
	for (int i = 1; i < MAX_EIGENVALUES * 2 && _Ψ[i].get() != NULL; i++)
		Ψ = Ψ + C[i] * _Ψ[i];
	} else {
		//define 2 test wave packets wave functions
		double σ = 0.2;
		double ωr = π;
		double ωi = π;
		ComplexFunctionPtr Ψ1(new ComplexFunction(X, [=](double x)->double{return gauss_cos(σ, 0.1*X, x, ωr); }, [=](double x)->double{return gauss_sin(σ, 0.1*X, x, ωi); }, X0, _precision));
		ComplexFunctionPtr Ψ2(new ComplexFunction(X, [=](double x)->double{return gauss_cos(σ, 0.5*X, x, ωr); }, [=](double x)->double{return gauss_sin(σ, 0.5*X, x, ωi); }, X0, _precision));

		Ψ = C[0] * Ψ1 + C[1] * Ψ2;//demo version
	}

	Ψ->Normalize();
	OperatorsTest(Ψ);
}

#ifndef _DEBUG
#define LAMBA_IMPLEMENTATION //use LAMBDA for SolveODE. WARNING: debug with lambda is very slow.
#endif

//define one of the standard potential(s) or none for the free particle
#define HARMONIC_OSCILATOR
//#define FINITE_POTENTIAL

//compute and output parameters
#define VALUE_FORMAT	"%lG"  //should be same digits as last step value (ex. 1E-5 -> %.5lf )


//universal constants
static const double ℏ = 1.054571726E-34; // reduced planck constant
static const double _eV=1.602176565E-19; // electron volt
static const double _me=9.10938291E-31; //electron mass
static const double _c = 299792458; //speed of light nm/s

static const double Hc = 2 / ℏ * _me / ℏ * _eV / 1E18; //this is schrodinger equation coeff 2m/ℏ^2 for electron mass and energy in eV, length in nm
/*
Solve time-independent schrodinger equation
p'(x)=-2m/h^2 ( En - V(x) ) Ψ(x)
Ψ'(x) = p(x)
where Ψ(x) is the wave function and V(x) is the time independent potential
*/



static double (*V)(double);

#define p F(0)
#define Ψ F(1)


#ifndef LAMBA_IMPLEMENTATION
/// Uncomment and replace SolveODE(odeFunc,... with SolveODENL(odefFuncNL,....
//p'(x)=-2m/h^2 ( En - V(x) ) Ψ(x)
double Energy(double x, FunctionSystem& F,void* udata)
{
	double En = *(double*)udata;
	return -Hc * m * (En - V(x)) * Ψ;
}

//Ψ'(x) = p(x)
double QMomentum(double x, FunctionSystem& F, void*)
{
	return p;
}
static FunctionSystem::FunctionDef odeFuncNL[2] = { Energy, QMomentum };
#endif

#ifdef HARMONIC_OSCILATOR //Harmonic oscilator
static double ω = 1E6; //ω in 1E+9 rad/s = 2*pi*freq (freq in GHz)
static const double _HarmonicC = (_me / _eV / 2); // 1/2 * me
static double _Harmonic_d0 = sqrt(ℏ *1E9 /*nm^2 / ps*/ / (m * _me * ω ));   //harmonic lengh = sqrt(ℏ / (m ω)) in nm
static double HARMONIC_POTENTIAL(double x){
	return (_HarmonicC * m * pow(ω*x, 2));
}
double HARMONIC_THEORETICAL_VAL(int odd, int i) {
	return	double(2 * odd + 4 * i + 1)*(ω*1.E9*ℏ / _eV/2); //(k+1/2)ℏω ,k>0
};
#endif

#ifdef FINITE_POTENTIAL //Finite Potential Step Bound States
static const double L = X/4; //potential well width in nm, the 1st step is at L/2
#define V0 10 //eV
double STEP_V(double x){ 
	return (x < L / 2) ? 0 : V0;
}
#define MIN_tan 1.E-4
#define MAX_tan 1.E4
double GetFiniteWellEn(int odd, int i,double EnMax) {
	static const int max_nodes = 2 * QuantumMechanic::MAX_EIGENVALUES;
	static double* even_nodes = NULL; //symetric solutions 
	static int even_count = 0;
	static double* odd_nodes; //anti-symetric solution
	static int odd_count = 0;
	static const double _C = 2 * ℏ / _me * ℏ / _eV * 1E18 / pow(L, 2); // (2ℏ^2) / (mL^2)  L in nm, energy in eV
	static const double _u0 = sqrt(EnMax / _C); // u0

	if (!even_nodes) { //solve eqns
		odd_nodes = new double[max_nodes];
		even_nodes = new double[max_nodes];
		Function oddf(_u0, [EnMax](double v){
			double tanv = tan(v);
			if (fabs(tanv) < MIN_tan) tanv = tanv > 0 ? MIN_tan : -MIN_tan;

			return sqrt(EnMax / _C - pow(v, 2)) + v / tanv;
		}, 0, 6);
		Function evenf(_u0, [EnMax](double v){
			double tanv = tan(v);
			if (fabs(tanv) > MAX_tan) tanv = tanv > 0 ? MAX_tan : -MAX_tan;
			return sqrt(EnMax / _C - pow(v, 2)) - v * tanv;
		}, 0, 6);
		AnalyzesResult odd_A{ odd_nodes, max_nodes };
		AnalyzesResult even_A{ even_nodes, max_nodes };


		evenf.Analyze(even_A);
		even_count = min(QuantumMechanic::MAX_EIGENVALUES, even_A.nodes);
		oddf.Analyze(odd_A);
		odd_count = min(QuantumMechanic::MAX_EIGENVALUES, odd_A.nodes);

		/*Function * funcToPlot[] = { &evenf, &oddf };
		const char* titles[] = { "solve finite well eigenvalue eqn",  NULL };
		const char* legends[] = { "evenf", "oddf",  NULL };
		Plot(sizeof(funcToPlot) / sizeof(funcToPlot[0]), funcToPlot, titles, "v", NULL, legends);*/
	}

	if (i >= (odd ? odd_count : even_count))
		return 0;

	return _C*pow((odd ? odd_nodes[i] : even_nodes[i]) /*v*/, 2);
}
#endif

double ZERO_POTENTIAL(double x){
	return  0;
}
void QuantumMechanic::BoundStates(const int n, bool SAVE_ALL_GRAPHS)
{
	int req_odd = -1, req_sol = -1;
	if (n >= 0) {
		req_odd = n % 2;
		req_sol = n / 2;
	}
//----------------- Config Start -------------------
	
#if defined HARMONIC_OSCILATOR
	auto THEORETICAL_VAL = HARMONIC_THEORETICAL_VAL;
	double EnMax = HARMONIC_POTENTIAL(X/4);
	V = HARMONIC_POTENTIAL;
#elif defined FINITE_POTENTIAL
	double EnMax = V0;
	auto THEORETICAL_VAL = [EnMax](int odd, int i)->double { return GetFiniteWellEn(odd,i,EnMax); };
	V = STEP_V;
#else
	double EnMax = 1;
	auto THEORETICAL_VAL = [](int odd, int i)->double { return 0; };
	V = ZERO_POTENTIAL;
#endif

//------------------ Config End --------------------------

//Computation params
	int const ENABLE_DBG = 0;  // 1 to enable [DBG] Output
	int const MIN_PRECISION = 1; //recommended 1 ,increase when expected compute values are less then 1E-MIN_PRECISION
	double const _step[] = { 1E-2, 1E-3, 1E-4
#ifndef _DEBUG
				,1E-5, 1E-6/*, 1E-7*/ // Release recommended { 1E-2, 1E-3, 1E-4, 1E-5 }, minim 2 values , add a 5-th step 1E-6 for high precision
#endif
	}; 
	int _MS = sizeof(_step) / sizeof(_step[0]);
	long long computationalTotal = 0; //total number of function called
//End Computation Params

//Start Computation
	unsigned int tick = os_ticks(), tickn;
	unsigned int first_tick = tick;
	double eigenValues[2][MAX_EIGENVALUES];
	int eigenValueCount[2];
	//preliminary step -- find local minimums
	if (n<0)
		printf("Start search for all eigenvalues between 0 and %.3lg\n", EnMax);
	else 
		printf("Start search for %d-th eigenvalue between 0 and %.3lg\n", n,EnMax);

	int computeTimes = (int)(EnMax / _step[0])+1;
	double* norm = new double[computeTimes];
	_precision = 1 + MIN_PRECISION - (int)round(log10(_step[0]));
	if (ENABLE_DBG) printf("\n[DBG][n]\tCandidate\n");
	for(int odd=((req_odd >=0) ? req_odd : 0);odd <= ((req_odd >=0) ? req_odd : 1);odd++)
	{
		double initialValues[2] = {
			(double)odd , // p[X0]
			1 - (double)odd	// y[X0]
		};
		computationalTotal += long long(8 * (X - X0) * Function::__Pow10(_precision)) * computeTimes; //8 is from 2 function and 4 calls from RK4
		#pragma omp parallel for
		for (int i = 0; i < computeTimes; i++)
		{
			double En = (i+1)*_step[0];
			FunctionSystem fs(2, X, initialValues, X0, -MIN_PRECISION);
			
#ifndef LAMBA_IMPLEMENTATION
			fs.SolveODENL(odeFuncNL, (void*)&En, FunctionSystem::Method::RK4, _precision);
#else
			std::function<double(double t, FunctionSystem& F)> odeFunc[2] = {
				[En](double x, FunctionSystem& F){ return -Hc * m * (En - V(x)) * Ψ; },
				[](double x, FunctionSystem& F){return p; }
			};
			fs.SolveODE(odeFunc, FunctionSystem::Method::RK4, _precision);
#endif
			norm[i] = fabs(*fs[1]);
		}

		eigenValueCount[odd] = 0;
		for (int i = 2; i < computeTimes; i++)
		{
			if (norm[i-2] > norm[i - 1] && norm[i - 1] < norm[i])
			{
				if (eigenValueCount[odd] == MAX_EIGENVALUES)
				{
					printf("Error: Not Enough Eigenvalue Buffer\n");
					break;
				}

				double minNormEn = i*_step[0];  //norm[i - 1] energy
				if (ENABLE_DBG) printf("[DBG][%d]\t" VALUE_FORMAT "\n", eigenValueCount[odd], minNormEn);
				eigenValues[odd][eigenValueCount[odd]++] = minNormEn;
			}
		}
	}
	delete norm;

	if (req_odd >= 0 && req_sol >= 0 && req_sol >= eigenValueCount[req_odd])
	{
		printf("Error requested solution not avaliable\n");
		return;
	}
	tickn = os_ticks();
	printf("-- Convergence analysis Step %d ready in %d ms -- \n", 0, tickn - tick);
	tick = tickn;
	//in-depth convergence steps -- find minimum with precision
	for(int s=1; s < _MS; s++)
	{
		_precision = MIN_PRECISION - (int)round(log10(_step[s]));
		computationalTotal += (req_sol >= 0 ? 1 : 2)/*odd and even*/ * long long(8 * (X - X0) * Function::__Pow10(_precision)) * computeTimes; //8 is from 2 function and 4 calls from RK4
		
		if (s == _MS - 1) {
			printf("\n[n]\tvalue\tteoretic value delta\n");
		}
		else if (ENABLE_DBG) {
			printf("\n[DBG][n]\tteoretic value delta\n");
		}
		
		for(int odd=((req_odd >=0) ? req_odd : 0);odd <= ((req_odd >=0) ? req_odd : 1);odd++)
		{
			double initialValues[2] = {
				(double)odd , // p[X0]
				1 - (double)odd	// y[X0]
			};

			int computeTimes = (int)((2* _step[s - 1]) / _step[s])-2;
			
			#pragma omp parallel for
			for(int i=((req_sol >=0) ? req_sol : 0);i<=((req_sol >=0) ? req_sol : eigenValueCount[odd]-1) ;i++)
			{
				double MinEn;
				double MinNorm=HUGE_VAL;
				
				#pragma omp parallel for
				for (int j = 0; j < computeTimes; j++)
				{
					double En = eigenValues[odd][i] - _step[s - 1] + (j+1) * _step[s];
					FunctionSystem fs(2, X, initialValues, X0, -MIN_PRECISION);
#ifndef LAMBA_IMPLEMENTATION
					fs.SolveODENL(odeFuncNL,(void*)&En, FunctionSystem::Method::RK4, _precision);
#else			
					std::function<double(double t, FunctionSystem& F)> odeFunc[2] = {
						[En](double x, FunctionSystem& F){ return -Hc * m * (En - V(x)) * Ψ; },
						[](double x, FunctionSystem& F){return p; }
					};
					fs.SolveODE(odeFunc, FunctionSystem::Method::RK4, _precision);
#endif

					double norm = fabs(*fs[1]);
					
					if (MinNorm > norm) {
							#pragma omp critical
							{
								if (MinNorm > norm) // compare again to avoid some other thread change the MinNorm in the mean time to a lower value
								{
									MinNorm = norm;
									MinEn = En;
								}
							}
					}
					
				}

				eigenValues[odd][i] = MinEn;

				if (ENABLE_DBG && s < _MS - 1)
				{
					#pragma omp critical
					{
						printf("[DBG][%d]\t" VALUE_FORMAT "\n", odd + 2 * i, MinEn - THEORETICAL_VAL(odd, i));
					}
				}

				if (s == _MS - 1)
				{
					#pragma omp critical
					{
						printf("[%d]\t" VALUE_FORMAT " eV\t" VALUE_FORMAT " eV\n", odd + 2 * i, MinEn, MinEn - THEORETICAL_VAL(odd, i));
					}
				}
			}
		}

		tickn = os_ticks();
		printf("-- Convergence analysis Step %d ready in %d ms -- \n", s, tickn - tick);
		tick = tickn;
	}
	
	int accuracy = max(0,_precision - 5);
	_precision = 5;
	if (req_odd >= 0 && req_sol >= 0) {
		computationalTotal += long long(8 * (X - X0) * Function::__Pow10(_precision + accuracy));
		if (eigenValueCount[req_odd] > req_sol)
		{
			double En = eigenValues[req_odd][req_sol];
			char title[256];
			snprintf(title, sizeof(title) - 1, "En=" VALUE_FORMAT " eV (odd %d, sol no. %d)", En, req_odd, req_sol);
			double initialValues[2] = {
				(double)req_odd, // p[X0]
				1 - (double)req_odd	// y[X0]
			};

			FunctionSystem fs(2, X, initialValues, X0, _precision);
			
#ifndef LAMBA_IMPLEMENTATION
			fs.SolveODENL(odeFuncNL, (void*)&En, FunctionSystem::Method::RK4, accuracy);
#else
			std::function<double(double t, FunctionSystem& F)> odeFunc[2] = {
				[En](double x, FunctionSystem& F){ return -Hc * m * (En - V(x)) * Ψ; },
				[](double x, FunctionSystem& F){return p; }
			};
			fs.SolveODE(odeFunc, FunctionSystem::Method::RK4, accuracy);
#endif
				
			printf("\n------ psi(x) ------ \n");

			AnalyzesResult a{};
			fs[1]->Analyze(a, stdout);
			double XMax = a._ConvergenceTime > 0 ? a._ConvergenceTime : X;
			printf("\n------ psi'(x) ----- \n");
			fs[0]->Analyze(a, stdout);

			fs[1]->ConvergeRight(XMax);
			long long N1 = fs[1]->CutRight(XMax);
			fs[0]->CutRight(XMax);

			Function PotentialFunc(XMax, V, X0);
			
			fs[1]->Normalize();
			fs[0]->Normalize();
			
			tickn = os_ticks();
			printf("\n---------------------\nCompute all ready int %u ms\nTOTAL COMPUTATION STEPS: %I64u\n---------------------------------\n", tickn - first_tick, computationalTotal);

			Function * funcToPlot[] = { fs[1], &PotentialFunc, fs[0] };
			const char* titles[] = { title, "", NULL };
			const char* legends[] = { "$\\psi(x)$", "V(x)", "$\\psi '(x)$", NULL };


			Plot(sizeof(funcToPlot) / sizeof(funcToPlot[0]), funcToPlot, titles, "x (nm)", NULL, legends);

			fs[1]->UndoCutRight(N1);
			_Ψ[req_odd + 2 * req_sol] = ComplexFunctionPtr(new ComplexFunction(fs.Release(1), NULL));
		}
	}
	else //save them all
	{
		if(SAVE_ALL_GRAPHS) {
			char outputDir[256];
			snprintf(outputDir, sizeof(outputDir) - 1, "Out_%08X",tick);
			CreateDirectoryA(outputDir, NULL);
			SetCurrentDirectoryA(outputDir);
		}
		
		for (int odd = 0; odd <= 1; odd++)
		{
			double initialValues[2] = {
				(double)odd, // p[X0]
				1 - (double)odd	// y[X0]
			};

			computationalTotal += long long(8 * (X - X0) * Function::__Pow10(_precision + accuracy) * eigenValueCount[odd]); //8 is from 2 function and 4 calls from RK4
			#pragma omp parallel for
			for (int i = 0; i <= eigenValueCount[odd] - 1; i++)
			{
				double En = eigenValues[odd][i];
				FunctionSystem fs(2, X, initialValues, X0,_precision);
				
#ifndef LAMBA_IMPLEMENTATION
				fs.SolveODENL(odeFuncNL, (void*)&En,FunctionSystem::Method::RK4);
#else
				std::function<double(double t, FunctionSystem& F)> odeFunc[2] = {
					[En](double x, FunctionSystem& F){ return -Hc * m * (En - V(x)) * Ψ; },
					[](double x, FunctionSystem& F){return p; }
				};
				fs.SolveODE(odeFunc,FunctionSystem::Method::RK4);
#endif
				AnalyzesResult a{};
				fs[1]->Analyze(a);
				double XMax = a._ConvergenceTime > 0 ? a._ConvergenceTime : X;
				fs[1]->ConvergeRight(XMax);
				fs[0]->ConvergeRight(XMax);
				
				fs[1]->Normalize();
				fs[0]->Normalize();
				if (SAVE_ALL_GRAPHS) {
					long long N1 = fs[1]->CutRight(XMax);
					fs[0]->CutRight(XMax);
					Function PotentialFunc(XMax, V, X0);
					Function * funcToPlot[] = { fs[1], &PotentialFunc, fs[0] };
					char title[256], title2[100];
					snprintf(title, sizeof(title) - 1, "Finite Well E%d.GIF", odd + 2 * i);
					snprintf(title2, sizeof(title2) - 1, "E%d = " VALUE_FORMAT " eV", odd + 2 * i, En);
					const char* titles[] = { title, title2, NULL };
					const char* legends[] = { "$\\psi(x)$", "V(x)", "$\\psi '(x)$", NULL };
					Plot(sizeof(funcToPlot) / sizeof(funcToPlot[0]), funcToPlot, titles, "x (nm)", NULL, legends, "GIF");
					fs[1]->UndoCutRight(N1);
				}
				_Ψ[odd+2*i] = ComplexFunctionPtr(new ComplexFunction(fs.Release(1), NULL));
			}
		}

		tickn = os_ticks();
		printf("\n---------------------\nCompute all ready int %u ms\nTOTAL COMPUTATION STEPS: %I64u\n---------------------------------\n", tickn - first_tick, computationalTotal);
	}
}
#undef p
#undef Ψ
//_ℏU is ℏ momentum operator -iℏ∂/∂x if x is in nm and p is multiple of eV/c
static const double _ℏU = ℏ / _eV/*eV*/ *_c*1E9 /*nm/s*/;

void QuantumMechanic::OperatorsTest(ComplexFunctionPtr& Ψ) {

	int tickn, tick = os_ticks();

	{//position operator x
		ComplexOperator<> Â{
			//real operator
			[](double x, ComplexFunction& f, long long s){
				return x*(f.real()[s]); 
			},
				[](double x, ComplexFunction& f, long long s){
				return x*(f.imag()[s]); 
			}
		};

		tick = os_ticks();

		ComplexFunctionPtr P = Â | Ψ;  // this is Â|Ψ>
		Complex expect_A_Ψ = Ψ | P;   // this is <Ψ|Â|Ψ>
		Complex expect_A2_Ψ = Ψ | (Â | P);   // this is <Ψ|Â^2|Ψ>

		tickn = os_ticks();
		printf("Compute in %d ms\n", tickn - tick);
		tick = tickn;

		char title2[128], title3[128];
		snprintf(title2, sizeof(title2) - 1, "$< \\psi |x| \\psi > =%.2lf + i(%.2lf) \\, nm$", expect_A_Ψ.real(), expect_A_Ψ.imag());
		snprintf(title3, sizeof(title3) - 1, "$\\Delta x=%.2lf + i($%.2lf)",
			expect_A2_Ψ.real() - (pow(expect_A_Ψ.real(), 2)), expect_A2_Ψ.imag() - (pow(expect_A_Ψ.imag(), 2)));


		const char* titles[] = { "Position Quantum Operator x", title2, title3, NULL };
		{
			Function*  funcToPlot[] = { Ψ.get()->Re(), P.get()->Re(), Ψ.get()->Im(), P.get()->Im() };
			const char* legends[] = { "$Re \\, \\psi$", "$Re \\, x | \\psi >$", "$Im \\, \\psi$", "$Im \\, x|\\psi >$", NULL };
			AnalyzesResult a{};
			for (int i = 0; legends[i]; i++)
			{
				printf("\n------ %s -------\n", legends[i]);
				funcToPlot[i]->Analyze(a, stdout);
			}
			Plot(sizeof(funcToPlot) / sizeof(funcToPlot[0]), funcToPlot,
				titles, "x nm", "", legends);
		}
	}

	{//Momentum Quantum Operator p
		
		ComplexOperator<> Â{
			//real operator
			[](double x, ComplexFunction& f, long long s){
				if (s < f.real().GetSize() - 1)
					return _ℏU*(f.imag()[s + 1] - f.imag()[s]) * f.real().GetUnitSteps(); //Re of momentum operator -iℏ∂/∂x
				return 0.;
			},
				[](double x, ComplexFunction& f, long long s){
				if (s < f.real().GetSize() - 1)
					return _ℏU * (-f.real()[s + 1] + f.real()[s]) * f.real().GetUnitSteps(); //Im of momentum operator -iℏ∂/∂x 
				return 0.;
			}
		};

		tick = os_ticks();

		ComplexFunctionPtr P = Â | Ψ;  // this is Â|Ψ>
		Complex expect_A_Ψ = Ψ | P;   // this is <Ψ|Â|Ψ>
		Complex expect_A2_Ψ = Ψ | (Â | P);   // this is <Ψ|Â^2|Ψ>

		tickn = os_ticks();
		printf("Compute in %d ms\n", tickn - tick);
		tick = tickn;

		char title2[128], title3[128];
		snprintf(title2, sizeof(title2) - 1, "$< \\psi |p| \\psi > =%.2lf + i(%.2lf) \\, nm$", expect_A_Ψ.real(), expect_A_Ψ.imag());
		snprintf(title3, sizeof(title3) - 1, "$\\Delta p=%.2lf + i($%.2lf)",
			expect_A2_Ψ.real() - (pow(expect_A_Ψ.real(), 2)), expect_A2_Ψ.imag() - (pow(expect_A_Ψ.imag(), 2)));


		const char* titles[] = { "Momentum Quantum Operator p", title2, title3, NULL };
		{
			Function*  funcToPlot[] = { P.get()->Re(), P.get()->Im() };
			const char* legends[] = { "$Re \\, p | \\psi > $", "$Im \\, p|\\psi > $", NULL };
			AnalyzesResult a{};
			for (int i = 0; legends[i]; i++)
			{
				printf("\n------ %s -------\n", legends[i]);
				funcToPlot[i]->Analyze(a, stdout);
			}
			Plot(sizeof(funcToPlot) / sizeof(funcToPlot[0]), funcToPlot,
				titles, "x nm", "", legends);
		}
	}

#ifdef HARMONIC_OSCILATOR
#define RAISE_TIMES 3
#define ANIHILATE_TIMES 2
	printf("\nLength of harmonic oscilator is: %lg nm\n\n",_Harmonic_d0);
	//Harmonic Oscillator - Raise Quantum Operator a+  1/sqrt(2) * (x/d0 + 1/d0 * ∂/∂x)
	{
		double c2 = (double)1 / sqrt(2);
		ComplexOperator<> Â{
			//real operator
			[c2](double x, ComplexFunction& f, long long s){
				if (s < f.real().GetSize() - 1)
					return c2*(x / _Harmonic_d0*f.real()[s] - _Harmonic_d0*(f.real()[s + 1] - f.real()[s]) * f.real().GetUnitSteps()); //Re of a+
				return 0.;
			},
				[c2](double x, ComplexFunction& f, long long s){
				if (s < f.imag().GetSize() - 1)
					return c2*(x / _Harmonic_d0*f.imag()[s] - _Harmonic_d0*(f.imag()[s + 1] - f.imag()[s]) * f.imag().GetUnitSteps()); //Im of a+
				return 0.;
			}
		};

		ComplexFunctionPtr P = ComplexFunctionPtr(new ComplexFunction(Ψ.get(),true)); //copy function
		for (int i = 1; i <= RAISE_TIMES; i++)
		{
			tick = os_ticks();
			P = Â | P; //P is Â^i|Ψ -- raise i times
			
			tickn = os_ticks();
			printf("Compute in %d ms\n", tickn - tick);
			tick = tickn;

			char title2[128];
			snprintf(title2, sizeof(title2) - 1, "raised %d time(s)", i);

			const char* titles[] = { "Harmonic Oscillator - Raise Operator a+", title2, NULL, NULL };
			{
				Function*  funcToPlot[] = { P.get()->Re(), P.get()->Im() };
				const char* legends[] = { "$Re \\, a+ | \\psi > $",  "$Im \\, a+|\\psi > $", NULL };
				AnalyzesResult a{};
				for (int i = 0; legends[i]; i++)
				{
					printf("\n------ %s -------\n", legends[i]);
					funcToPlot[i]->Analyze(a, stdout);
				}
				Plot(sizeof(funcToPlot) / sizeof(funcToPlot[0]), funcToPlot,
					titles, "x nm", "", legends);
			}
		}
	}

	//Harmonic Oscillator - Anihilator Quantum Operator a  1/sqrt(2) * (x/d0 + 1/d0 * ∂/∂x)
	{
		double c2 = (double)1 / sqrt(2);
		ComplexOperator<> Â{
			//real operator
			[c2](double x, ComplexFunction& f, long long s){
				if (s < f.real().GetSize() - 1)
					return c2*(x / _Harmonic_d0*f.real()[s] + _Harmonic_d0*(f.real()[s + 1] - f.real()[s]) * f.real().GetUnitSteps()); //Re of a
				return 0.;
			},
				[c2](double x, ComplexFunction& f, long long s){
				if (s < f.imag().GetSize() - 1)
					return c2*(x / _Harmonic_d0*f.imag()[s] + _Harmonic_d0*(f.imag()[s + 1] - f.imag()[s]) * f.imag().GetUnitSteps()); //Im of a
				return 0.;
			}
		};

		ComplexFunctionPtr P = ComplexFunctionPtr(new ComplexFunction(Ψ.get(), true)); //copy function
		for (int i = 1; i <= ANIHILATE_TIMES; i++)
		{
			tick = os_ticks();
			P = Â | P; //P is Â^i|Ψ -- anihilate i times

			tickn = os_ticks();
			printf("Compute in %d ms\n", tickn - tick);
			tick = tickn;

			char title2[128];
			snprintf(title2, sizeof(title2) - 1, "anihilated %d time(s)", i);

			const char* titles[] = { "Harmonic Oscillator - Anihilator Operator a", title2, NULL, NULL };
			{
				Function*  funcToPlot[] = { P.get()->Re(), P.get()->Im() };
				const char* legends[] = { "$Re \\, a | \\psi > $", "$Im \\, a|\\psi > $", NULL };
				AnalyzesResult a{};
				for (int i = 0; legends[i]; i++)
				{
					printf("\n------ %s -------\n", legends[i]);
					funcToPlot[i]->Analyze(a, stdout);
				}
				Plot(sizeof(funcToPlot) / sizeof(funcToPlot[0]), funcToPlot,
					titles, "x nm", "", legends);
			}
		}
	}
#endif

}


