#pragma once
#ifndef __FunctionalMath__943hjf3489yr43__
#define __FunctionalMath__943hjf3489yr43__

#include "target.h"
#include <functional>
#include <memory>
#include <complex>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif


#define double_FORMAT  "%lg"

namespace FunctionalMath {
	class Function;
	struct AnalyzesResult{
		//input
		//solve f(t)=0 and save t into nodesTime (max nodesTimeSize)
		double* nodesTime;
		int nodesTimeSize;
		//output
		double _MinVal, _min, _max, _MaxVal;
		int nodes, localMins, localMaxs;
		double _ConvergenceTime, _ConvergenceDelta;
	};

	struct Operator {
		std::function < double(double, Function&, long long) > op;
	};

	template<typename type = ComplexFunction>
	struct ComplexOperator {
		std::function < double(double, type&, long long) > opRe;
		std::function < double(double, type&, long long) > opIm;
	};
	using Complex = std::complex < double > ;
	class Function
	{
		friend class FunctionSystem;
		friend class ComplexFunction;
		friend class PDEData;
		Function() {} //private constructor used by internal functions

		double _t0; //interval start
		double _T; //interval stop
		double _UnitSteps;//compute steps per unit ;this is 10^_precision
		long long _N; // data points
		long long _lastIndex;

		//single args
		double _currentRealValue;
		double* data;

		//multi-args private members
		Function** fdata;
		int _level; //current function level
		void Construct(int level, int masterLevel, double* T, double *t0, int* precision,
			std::function<double(double* X)>* gen = NULL,
			double* X = NULL);

		static void ComputePartialDerivative(int k, Function* df, Function* f, Function* f1 = NULL, double unitSteps = 0);
		static double null_value_ref;
		static Function null_function_ref;

		

	public:
		class PDEData;
		typedef double(*PDEEqn)(double* X, PDEData& f, void*);
		class PDEData {
			std::random_device rd;
			std::mt19937 rgen;
			std::uniform_real_distribution<> rdistr;

			friend class Function;
			//global data
			Function& f;
			double* X; //current vector in real values
			long long* s;//current vector in discrete index values
			double* UnitSteps;//
			double valMin, valMax;
			double valSteps;
			enum __state : unsigned int {
				initial=0,
				select,
				compute
			}state;
			int n;//number of variables
			PDEEqn eqn;
			void* userdata;
			
			//adaptive params and stats
			double score; //aprox score -- 0 is the best
			double maxDelta; //variation absolute max
			double rDelta;//variation = [rDelta-1,2-rDelta]*maxDelta  
			double retriesCount;
			unsigned int maxRetries;
			unsigned int selectSteps;
			unsigned int computeSteps;

			//current step shared data
			double threshold;
			int unkMaxSize;
			double** unk1; //current unknowns references
			double* unk2; //current unknowns best values -- select phase
			double* unk3; //current unknowns best values -- compute phase
			int unkGlobalIndex; //current unknowns numbers (including duplicates)

			//current step per omp thread data
			int* unkIndex;//curent unknown indexes of all thread
			long long* omps;    //curent step vector(s) for all threads
			double* unkd; //curent unknowns values of all thread (no duplicates)

			//current step per omp thread functions
			void InitOMPData();
			void ResetOMPData();
			double SetDuplicateOMPUnk(int i);
			double* GetOMPUnk();
			void SetOMPUnk( double val);
			int GetOMPIndex();
			long long* GetOMPStepVector();
			bool boundary;	
			

			//global functions
			PDEData(Function&func) :f(func) { rgen.seed(rd()); } //internal class only cannot be instanciated outside Function
			~PDEData() {
				delete[] X;
				delete[] s;
				delete[] UnitSteps;
				delete[] unk1;
				delete[] unk2;
				delete[] unk3;
				delete[] unkd;
				delete[] omps;
				delete[] unkIndex;
			}
			double getRandomVal(double& fval, double* fref = NULL, double unitSteps = 0, double* fref2 = NULL, double unitSteps2 = 0);
			
		public:
			//current value of function ∂^2f/∂xi∂xj during solving PDE
			//for ∂f/∂xi use j=-1
			//for ∂^2f/∂xi^2 use j==i
			double operator()(int i, int j = -1);

			
		};
		

		Function** SolvePDENL(int& maxSolutions, PDEEqn, void* userdata, double delta, unsigned int selectSteps, unsigned int computeSteps, int precision, double rDelta = 0,double valMin = 0, double valMax = 0);
	private:
		bool SolvePDENL_MOLRecursive(int l, PDEData& pdeData);

	
	public:
		static long long maxN;
		static int maxNOutputLines;
		static bool DBG_PDE_SOLVE;


		static inline bool IsNull(double& val) { return &val == &null_value_ref; }
		static inline bool IsNull(Function& val) { return &val == &null_function_ref; }
		static double __Pow10(int n);

		//copy constructor(s)
		Function(const Function *f, bool bCopyData = false);
		Function(const Function& f) : Function(&f, true) {}
		inline Function& operator=(const Function& f) { Function(&f, true); }

		//single arg constructors
		
		Function(double T, double Y0 = 0, double t0 = 0, int precision = 4);
		Function(double T, std::function<double(double)> gen, double t0 = 0, int precision = 4);

		//multi-arg F(X) X={X0, ...,Xn-1}
		//n is the number of arguments
		//XMax is {X0Max, ...,Xn-1Max}
		//XMin is {X0Min, ...,Xn-1Min}
		//precision is {P0,P1,..Pn-1} is precision in each variable
		Function(int n, double* XMax, double *XMin, int* precision) {
			Construct(n, n, XMax, XMin, precision);
		}

		//multi-arg F(X) X={X0, ...,Xn-1}
		//n is the number of arguments
		//XMax is {X0Max, ...,Xn-1Max}
		//XMin is {X0Min, ...,Xn-1Min}
		//gen generates f(X) where X={X0, ...,Xn-1}
		//precision is {P0,P1,..Pn-1} is precision in each variable
		Function(int n, double* XMax, double *XMin, std::function<double(double* X)> gen, int* precision){
			Construct(n, n, XMax, XMin, precision,  &gen);
		}
		virtual ~Function(void);


		//fast calls
		inline double GetT0() { return _t0; }
		inline double GetT() { return _T; }
		inline long long GetSize() { return _N; }
		inline double GetUnitSteps() { return _UnitSteps; }
		inline unsigned int GetPrecision() { return (unsigned int)round(log10(_UnitSteps)); }
		inline int GetLevel() { return _level; }
		inline operator double*() { return data; }
		inline operator double() { return _currentRealValue; }
		inline void operator=(double val) { _currentRealValue = val; }
		inline void operator+=(double val) { _currentRealValue += val; }

		inline void SetAt(long long i, double val) { data[i] = val; }
		inline void ResetLastIndex() { _lastIndex = 0; _currentRealValue = data[0]; }
		inline double PopLastValue() { _currentRealValue = data[_lastIndex++]; return _currentRealValue; }

		inline void SaveCurrentValue() {
			if (_N - 1 > _lastIndex)
				data[++_lastIndex] = _currentRealValue;
		}

		inline void ResetData() {
			if (data)
				memset(data, 0x00, _N * sizeof(double));
		}
		

		//fast discrete operators -- no checking
		inline double& operator[](long long n) { return data[n]; }
		inline Function& operator[](int n) { return *fdata[n]; }

		inline double& operator[](long long* s){
			if (data) {
				if (*s < _N)
					return data[*s];
				else 
					return null_value_ref;
			}
			if (*s < _N)
				return (*fdata[*s])[s + 1];
			else
				return null_value_ref;
		}

		
		inline Function& operator()(long long* s){
			if (data) {
				return *this;
			}
			return (*fdata[*s])(s + 1);
		}

		//return sub-function f({xk+1,...,xn-1}) 
		//X={X0, ..., Xk}, k < n-1
		inline Function& operator()(double* X){
			if (data) {
				return *this;
			}

			long long n = (int)((*X - _t0)*_UnitSteps);
			if (n < 0) return null_function_ref;
			else if (n >= _N) return null_function_ref;
			else return (*fdata[n])(X + 1);
		}
		
		//return value f(X) where X={X0, ...,Xn-1}
		inline double& operator[](double* X) {
			if (data) return operator[](*X);

			long long n = (int)((*X- _t0)*_UnitSteps);
			if (n < 0) return null_value_ref;
			else if (n >= _N) return null_value_ref;
			else return (*fdata[n])[X + 1];
		}

		inline double& operator[](double t) {
			if (!data) return null_value_ref;

			long long n = (int)((t - _t0)*_UnitSteps);
			if (n < 0) return null_value_ref;
			else if (n >= _N) return null_value_ref;
			else return data[n];
		}

		//set all data to val
		void Reset(double val = 0);

		//converge constant between TMax and T with same value as f(TMax)
		void ConvergeRight(double TMax);
		
		//use to set T to TMax
		long long CutRight(double TMax);

		//inverse of CutRight _N to old one returned by CutRight
		//use the value returned by CutRight as parameter oldN
		void UndoCutRight(long long oldN);

		//generate a function with real part this Function and imaginary part 0
		operator ComplexFunction*();
		
		//calculate fourier
		ComplexFunction* FourierTransform(double MaxF, int precision = 2);

		//create a function I(t) where I(t) is integral over[t0,t] of this function..does not work on multiple args
		Function * CalculateIntegral();

		//create n-th order derivate function of this function..does not work on multiple args 
		Function * CalculateDerivate(unsigned int n = 1);

		//create n-th order partial derivate of Xk of f(X) where X={X0, ...,xk,...,Xn-1} 0 <= k < n
		Function * CalculatePartialDerivate(unsigned int k,unsigned int n = 1);

		//integrate on the entire interval
		inline double Integrate() { return Integrate(_t0, _T); }

		//integrate on the [t1,t2] interval. .does not work on multiple args
		double Integrate(double t1, double t2);

		//volume integration of multi args F({x0,..xn-1}) on {[a0,b0],...,[an-1,bn-1]}
		//N is the number of points in Monte-Carlo method
		double Integrate(double* a, double* b, int N= 1E6);

		//Analyze min,max,nodes -- print result on stdout
		void Analyze(AnalyzesResult& result, FILE* file = NULL, bool bReset = true);

		//Dump all data one per row to console
		void Dump(FILE* file = stdout);

		//normalize the function -- divide by sqrt((f*f)->Integral())
		void Normalize();
		
		

		//1st parameter of op is current double of t
		//2nd parameter of op is this function double, use only doubles of function before or including step
		//3rd paramter of op is current calculations step
		Function* ApplyOperator(const Operator& O);
		ComplexFunction* ApplyOperator(const ComplexOperator<Function>& O);

		
		typedef double(*OperatorDef)(double, Function&, long long);
		Function* ApplyOperatorNL(OperatorDef op);
		ComplexFunction* ApplyComplexOperatorNL(OperatorDef op_Re, OperatorDef op_Im);

		//inner and scalar products -- do not use this just use the FunctionPtr operators below
		Function* operator*(const Function*);
		Function* operator*(double);

		Function* operator/(const Function*);
		Function* operator/(double);

		Function* operator+(const Function*);
		Function* operator+(double);

		Function* operator-();
		Function* operator-(const Function*);
		Function* operator-(double);
		Function* SubstractFrom(double x);
	};
	
	
	using FunctionPtr = std::unique_ptr <Function>;
	inline FunctionPtr operator*(const FunctionPtr& f1, const FunctionPtr& f2) { return FunctionPtr((*f1) * f2.get()); }
	inline FunctionPtr operator*(double x, const FunctionPtr& f1) { return FunctionPtr(*f1*x); }
	inline FunctionPtr operator*(const FunctionPtr& f1, double x) { return FunctionPtr(*f1*x); }

	inline FunctionPtr operator/(const FunctionPtr& f1, const FunctionPtr& f2) { return FunctionPtr((*f1) / f2.get()); }
	inline FunctionPtr operator/(const FunctionPtr& f1, double x) { return FunctionPtr(*f1/x); }
	
	inline FunctionPtr operator+(const FunctionPtr& f1, const FunctionPtr& f2) { return FunctionPtr((*f1) + f2.get()); }
	inline FunctionPtr operator+(double x, const FunctionPtr& f1) { return FunctionPtr(*f1+x); }
	inline FunctionPtr operator+(const FunctionPtr& f1, double x) { return FunctionPtr(*f1+x); }
	
	inline FunctionPtr operator-(const FunctionPtr& f1) { return FunctionPtr(-(*f1)); }
	inline FunctionPtr operator-(const FunctionPtr& f1, const FunctionPtr& f2) { return FunctionPtr((*f1) - f2.get()); }
	inline FunctionPtr operator-(const double x, const FunctionPtr& f1) { return FunctionPtr(f1->SubstractFrom(x));}
	inline FunctionPtr operator-(const FunctionPtr& f1, double x) { return FunctionPtr(*f1 + x); }

	//this is O|f1> in Dirac notation
	inline FunctionPtr operator|(const Operator& O, FunctionPtr& f1) { return FunctionPtr((*f1).ApplyOperator(O)); }


	class ComplexFunction
	{
		Function* _Re;
		Function* _Im;
	public:
		inline Function* Re() { return _Re; }
		inline Function* Im() { return _Im; }
		inline Function& real() { return *_Re; }
		inline Function& imag() { return *_Im; }

		ComplexFunction(const ComplexFunction *f, bool bCopyData = false);
		ComplexFunction(const ComplexFunction& f) : ComplexFunction(&f, true) {}
		inline ComplexFunction& operator=(const ComplexFunction& f) { ComplexFunction(&f, true); }

		ComplexFunction(double T, double Re0 = 0, double Im0 = 0, double t0 = 0, int precision = 4);
		ComplexFunction(Function* f);
		ComplexFunction(double T, std::function<double(double)> genRe, std::function<double(double)> genIm,
						double t0 = 0, int precision = 4);
		ComplexFunction(Function* real, Function* imag);
		~ComplexFunction();

		Function* InverseFourierTransform(double T, double t0 = 0, int precision = 4);

		Function* Module();
		Function* Arg();
		
		//normalize the function -- divide by sqrt(<f|f>->Integral()) so the new function will have <f|f>->Integral = 1
		void Normalize();
		
		
		//1st parameter of op is current double of t
		//2nd parameter of op is this function double, use only doubles of function before or including step
		//3rd paramter of op is current calculations step
		//use NL version to reduce process time
		ComplexFunction* ApplyOperator(const ComplexOperator<>& O);

		typedef double(*OperatorDef)(double, ComplexFunction&, long long);
		ComplexFunction* ApplyComplexOperatorNL(OperatorDef op_Re, OperatorDef op_Im);

		
		ComplexFunction* operator*(const Function*);
		ComplexFunction* operator*(const ComplexFunction*);
		ComplexFunction* operator*(double x);
		ComplexFunction* operator*(const Complex& x);

		ComplexFunction* operator/(const Function*);
		ComplexFunction* operator/(const ComplexFunction*);
		ComplexFunction* operator/(double x);
		ComplexFunction* operator/(const Complex& x);

		ComplexFunction* operator+(const Function*);
		ComplexFunction* operator+(const ComplexFunction*);
		ComplexFunction* operator+(double x);
		ComplexFunction* operator+(const Complex& x);

		ComplexFunction* operator-();
		ComplexFunction* operator-(const Function*);
		ComplexFunction* operator-(const ComplexFunction*);
		ComplexFunction* operator-(double);
		ComplexFunction* operator-(const Complex&);
		static ComplexFunction* ScalarMinusFunction(double x, ComplexFunction& f1);
		static ComplexFunction* ScalarMinusFunction(const Complex& x, ComplexFunction& f1);
		static ComplexFunction* ScalarMinusFunction(const Function& x, ComplexFunction& f1);

		//scalar product (inner product)
		Complex operator|(const ComplexFunction*);
		ComplexFunction* operator~();
	};

	using ComplexFunctionPtr = std::unique_ptr <ComplexFunction>;
	
	inline ComplexFunctionPtr operator*(const ComplexFunctionPtr& f1, const FunctionPtr& f2) { return ComplexFunctionPtr((*f1) * f2.get()); }
	inline ComplexFunctionPtr operator*(const ComplexFunctionPtr& f1, const ComplexFunctionPtr& f2) { return ComplexFunctionPtr((*f1) * f2.get()); }

	inline ComplexFunctionPtr operator*(const FunctionPtr& f1, const ComplexFunctionPtr& f2) { return ComplexFunctionPtr((*f2) * f1.get()); }
	inline ComplexFunctionPtr operator*(double x, const ComplexFunctionPtr& f1) { return ComplexFunctionPtr(*f1*x); }
	inline ComplexFunctionPtr operator*(const Complex& x, const ComplexFunctionPtr& f1) { return ComplexFunctionPtr(*f1*x); }
	inline ComplexFunctionPtr operator*(const ComplexFunctionPtr& f1, double x) { return ComplexFunctionPtr(*f1*x); }
	inline ComplexFunctionPtr operator*(const ComplexFunctionPtr& f1, const Complex& x) { return ComplexFunctionPtr(*f1*x); }

	inline ComplexFunctionPtr operator/(const ComplexFunctionPtr& f1, const FunctionPtr& f2) { return ComplexFunctionPtr((*f1) / f2.get()); }
	inline ComplexFunctionPtr operator/(const ComplexFunctionPtr& f1, const ComplexFunctionPtr& f2) { return ComplexFunctionPtr((*f1) / f2.get()); }
	inline ComplexFunctionPtr operator/(const ComplexFunctionPtr& f1, double x) { return ComplexFunctionPtr(*f1/x); }
	inline ComplexFunctionPtr operator/(const ComplexFunctionPtr& f1, const Complex& x) { return ComplexFunctionPtr(*f1/x); }

	inline ComplexFunctionPtr operator+(const ComplexFunctionPtr& f1, const FunctionPtr& f2) { return ComplexFunctionPtr((*f1) + f2.get()); }
	inline ComplexFunctionPtr operator+(const FunctionPtr& f1, const ComplexFunctionPtr& f2) { return ComplexFunctionPtr((*f2) + f1.get()); }
	inline ComplexFunctionPtr operator+(const ComplexFunctionPtr& f1, const ComplexFunctionPtr& f2) { return ComplexFunctionPtr((*f1) + f2.get()); }
	inline ComplexFunctionPtr operator+(double x, const ComplexFunctionPtr& f1) { return ComplexFunctionPtr(*f1+x); }
	inline ComplexFunctionPtr operator+(const ComplexFunctionPtr& f1, double x) { return ComplexFunctionPtr(*f1+x); }
	inline ComplexFunctionPtr operator+(const Complex& x, const ComplexFunctionPtr& f1) { return ComplexFunctionPtr(*f1 + x); }
	inline ComplexFunctionPtr operator+(const ComplexFunctionPtr& f1, const Complex& x) { return ComplexFunctionPtr(*f1 + x); }

	inline ComplexFunctionPtr operator-(const ComplexFunctionPtr& f1, const FunctionPtr& f2) { return ComplexFunctionPtr((*f1) - f2.get()); }
	inline ComplexFunctionPtr operator-(const FunctionPtr& f1, const ComplexFunctionPtr& f2) { return ComplexFunctionPtr(ComplexFunction::ScalarMinusFunction(*f1,*f2)); }
	inline ComplexFunctionPtr operator-(const ComplexFunctionPtr& f1) { return ComplexFunctionPtr(-(*f1)); }
	inline ComplexFunctionPtr operator-(const ComplexFunctionPtr& f1, const ComplexFunctionPtr& f2) { return ComplexFunctionPtr((*f1) - f2.get()); }
	inline ComplexFunctionPtr operator-(double x, const ComplexFunctionPtr& f1) { return ComplexFunctionPtr(ComplexFunction::ScalarMinusFunction(x, *f1));  }
	inline ComplexFunctionPtr operator-(const ComplexFunctionPtr& f1, double x) { return ComplexFunctionPtr(*f1 + x); }
	inline ComplexFunctionPtr operator-(const Complex& x, const ComplexFunctionPtr& f1) { return ComplexFunctionPtr(ComplexFunction::ScalarMinusFunction(x, *f1)); }
	inline ComplexFunctionPtr operator-(const ComplexFunctionPtr& f1, const Complex& x) { return ComplexFunctionPtr(*f1 + x); }
	
	inline ComplexFunctionPtr operator~(const ComplexFunctionPtr& f1) { return ComplexFunctionPtr(~(*f1)); }
	//Dirac notation f| is the bra <f| and |f is the ket |f>

	//this is inner product <f1|f2> in Dirac notation
	inline Complex operator|(const ComplexFunctionPtr& f1, const ComplexFunctionPtr& f2)     { return (*f1)|(f2.get()); }
	//this is O|f1> in Dirac notation
	inline ComplexFunctionPtr operator|(const ComplexOperator<ComplexFunction>& O, const ComplexFunctionPtr& f1)  { return ComplexFunctionPtr((*f1).ApplyOperator(O)); }
	//this is O|f1> in Dirac notation
	inline ComplexFunctionPtr operator|(const ComplexOperator<Function>& O, const FunctionPtr& f1) { return ComplexFunctionPtr((*f1).ApplyOperator(O)); }

	class FunctionSystem
	{
	public:
		class ODEAnalyzesResult{
			friend class FunctionSystem;
		public:
			double _MinVal, _min, _max, _MaxVal;
			int nodes;
			double _ConvergenceTime, _ConvergenceDelta;
			//compute internal data
		private:
			double _ConvergenceStart, _ConvergenceEnd, _Convergencedouble;
			double _prevdouble;
		};
		typedef double(*FunctionDef)(double t, FunctionSystem& F,void* data);
		typedef double(*FunctionMDef)(double* X, FunctionSystem& F, void* data);

	private:
		int _n;
		Function **_Y;
		void Init(ODEAnalyzesResult*res, double t0, double initialVal);

		void Compute(ODEAnalyzesResult*res, double t, double val);

	public:
		/*create a F0,F2,...Fn-1 functions system on interval [t0,T]*/
		FunctionSystem(int n, double T, double* initialdoubles, double t0 = 0, int precision = 4);
		
		//create a F0,F2,...Fn-1 functions system of m variables
		//n is the number of arguments
		//XMax is {X0Max, ...,Xm-1Max}
		//XMin is {X0Min, ...,Xm-1Min}
		//precision is {P0,P1,..Pm-1} is precision in each variable
		//gen(X,k) is used to specifiy initial condition for Fk
		//gen should clearly generate all Fk(XMin0,XMin1,...,xi,...,XMin[m-1]) where k=0..n-1 and i=0..m
		//gen should return HUGE_VAL for unknown data
		FunctionSystem(int n, int m, double* XMax, double *XMin, std::function<double(double* X,int k)> gen, int* precision);

		~FunctionSystem(void);

		inline Function* operator[](int n) { return _Y[n]; }

		//release n-th function from the function system operator[n] will return NULL from now on and the function will not be deleted.
		//return the n-th function pointer
		inline Function* Release(int n) { Function* f = _Y[n]; _Y[n] = NULL; return f; }

		// current value of function Fk during solving ODE
		inline double operator()(int k) { 
			return (*_Y[k]);
		}

		//compute function methods
		enum Method{
			RK4
		};



		
		/*Solve ODE eqn system Fk'(t)=f[k](t,F0(t),F2(t),..Fn-1(t)) for k=0..n-1
		- f[k] is defined as a set of k=0..n-1 FunctionDef(s)
		- use F(k) for Fk(t) inside each f[k] function implementation
		- userdata will be passed on to each functionDef call
		*/
		bool SolveODENL(FunctionDef* f, void* userdata = NULL, Method method = RK4, unsigned int accuracy = 0, ODEAnalyzesResult* res = NULL);

		/*Solve ODE eqn system Fk'(t)=f[k](t,F0(t),F2(t),..Fn-1(t)) for k=0..n-1
		- f[k] is defined as a set of k=0..n-1 FunctionDef(s)
		- use F(k) for Fk(t) inside each f[k] function implementation
		WARNING: debug with lambda is very slow
		*/
		bool SolveODE(std::function<double(double t, FunctionSystem& F)>* f, Method method = RK4, unsigned int accuracy = 0, ODEAnalyzesResult* res = NULL);
		
		/*Eval function based on current function system
		- use F(k) for Fk(t)
		- NOTE: for faster calculation use NL version*/
		Function* CalculateFunction(std::function<double(double, FunctionSystem& F)> gen);
		Function * CalculateFunctionNL(FunctionDef gen,void*userdata=NULL);
	};



};
#endif