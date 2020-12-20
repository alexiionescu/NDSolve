#include "Function.h"
using namespace FunctionalMath;
#ifndef π
#define π 3.14159265358979323846
#endif

long long Function::maxN = (long long)(1.E+9); // global max memory (multiple of 8) for all Function instances. 
double Function::null_value_ref = 0;
Function Function::null_function_ref{};
bool Function::DBG_PDE_SOLVE = false;

double Function::__Pow10(int n){
	static double pow10[11] = {
		1, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6, 1E7, 1E8, 1E9, 1E10 };
	static double pow_10[11] = {
		1, 1E-1, 1E-2, 1E-3, 1E-4, 1E-5, 1E-6, 1E-7, 1E-8, 1E-9, 1E-10 };
	double x;
	if (n < 0)
	{
		x = pow_10[n > -10 ? -n : 10];
		while (n++ < -10) x /= 10;
	}
	else
	{
		x = pow10[n < 10 ? n : 10];
		while (n-- > 10) x *= 10;
	}

	return x;
}

Function::Function(double T, double Y0, double t0, int precision)
{

	_UnitSteps = __Pow10(precision);
	long long totalN = long long((T - t0)*_UnitSteps);
	if (totalN > maxN) {
		int p = (int)ceil(log10(totalN / maxN));
		precision -= p;
		_UnitSteps = __Pow10(precision);
	}
	
	_t0 = t0;
	_T = T;
	_N = (long long)((T - t0)*_UnitSteps + 1);
	if (_N < 2) _N = 2; // at least initial double and last double are kept
#pragma omp critical
	{
		maxN -= sizeof(Function) / sizeof(double);
		maxN -= _N;
	}
	if (maxN < 0) throw "Out of Memory Exception";


	data = new double[_N];
	_level = 1;
	_currentRealValue = *data = Y0;
	_lastIndex = 0;
	fdata = NULL;
}

Function::Function(const Function *f, bool bCopyData)
{
	_T = f->_T;
	_t0 = f->_t0;
	_UnitSteps = f->_UnitSteps;
	_N = f->_N;
#pragma omp critical
	{
		maxN -= sizeof(Function) / sizeof(double);
		maxN -= _N;
	}
	if (maxN < 0) throw "Out of Memory Exception";
	_level = f->_level;
	_lastIndex = 0;

	if (f->data) {
		data = new double[_N];
		
		fdata = NULL;
		_currentRealValue = *data = f->data[0];
		
		if (bCopyData)
			memcpy(data + 1, f->data + 1, (_N - 1)*sizeof(double));
	}
	else
	{
		data = NULL;
		fdata = new Function*[_N];
		
		for (long long i = 0; i < _N; i++) {
			fdata[i] = new Function(f->fdata[i],bCopyData);
		}
	}
	
}


Function::Function(double T, std::function<double(double)> gen, double t0, int precision)
{
	_UnitSteps = __Pow10(precision);
	long long totalN = (long long)((T - t0)*_UnitSteps);
	if (totalN > maxN) {
		int p = (int)ceil(log10(totalN / maxN));
		precision -= p;
		_UnitSteps = __Pow10(precision);
	}
	_t0 = t0;
	_T = T;
	_N = (long long)((T - t0)*_UnitSteps + 1);
	if (_N < 2) _N = 2; // at least initial  and last value are kept
#pragma omp critical
	{
		maxN -= sizeof(Function) / sizeof(double);
		maxN -= _N;
	}
	if (maxN < 0) throw "Out of Memory Exception";

	data = new double[_N];

	for (int j = 0; j < _N; j++) {
		data[j] = gen(_t0+j/_UnitSteps);
	}
	_lastIndex = 0;
	_level = 1;
	_currentRealValue = data[0];
	fdata = NULL;
}

//multi-arg F(X) X={X0, ...,Xn-1}
//MasterLevel is the number of arguments
//level is current built sub-function level. starts equalwith MasterLevel -> 1
//XMax is {X0Max, ...,Xn-1Max}
//XMin is {X0Min, ...,Xn-1Min}
//Y0 generates F(0,0,...,Xn-1)
//gen generates entire F(X) when level reaches 1
void Function::Construct(int level, int n, double* XMax, double *XMin, int* precision,
					std::function<double(double* X)>* gen,
					double* X)
{
	_level = level;
	if (n == level) {
		long long totalN = 1;
		for (int i = 0; i < n; i++) {
			totalN *= (long long)((XMax[i] - XMin[i])*__Pow10(precision[i]));
			if (i < n - 1) totalN *= sizeof(Function) / sizeof(double);
		}

		if (totalN > maxN) {
			int p = (int)max(1,ceil(log10(totalN / maxN)/level));
			for (int i = 0; i < n; i++)
				precision[i] -= p;
		}

	}

	_UnitSteps = __Pow10(precision[n - level]);
	_t0 = XMin[n - level];
	_T = XMax[n - level];
	_N = (long long)((_T - _t0)*_UnitSteps + 1);
	if (_N < 2) _N = 2; // at least initial  and last value are kept
#pragma omp critical
	{
		maxN -= sizeof(Function) / sizeof(double);
		maxN -= _N;
	}
	if (maxN < 0) throw "Out of Memory Exception";
	_lastIndex = 0;
	
	if (level == 1)	{ //Xn-1 level
		data = new double[_N];
		
		fdata = NULL;
		if (gen) {
			for (int j = 0; j < _N; j++) {
				X[n - 1] =  _t0 + j / _UnitSteps;
				data[j] = (*gen)(X);
			}
			_currentRealValue = *data;
		}
	}
	else { //function values

		data = NULL;
		fdata = new Function*[_N];
		
		if (n == level && gen)
			X = new double[level];

		for (long long j = 0; j < _N; j++) {
			fdata[j] = new Function();
			if (gen) X[n - level] = _t0 + j / _UnitSteps;
			fdata[j]->Construct(level - 1, n, XMax, XMin, precision, gen,  X);
		}

		if (n == level && gen)
			delete[] X;
	}
}

Function::~Function(void)
{
#pragma omp critical
	{
		maxN += sizeof(Function) / sizeof(double);
		maxN += _N;
	}

	if (data) {
		delete[] data;
	}

	if (fdata) {
		for (long long i = 0; i < _N; i++)
			delete fdata[i];
		delete[] fdata;
	}
	
}

//calculate fourier
ComplexFunction* Function::FourierTransform(double MaxF, int precision){
	if (!data) return NULL;
	ComplexFunction* F = new ComplexFunction(MaxF, 0, 0, 0, precision);
	long long computesteps = F->Re()->GetSize();
	const double _2pi = 2 * π;

#pragma omp parallel for
	for (long long i = 0; i < computesteps; i++)
	{
		double _re = 0, _im = 0;
		double _2pi_f = _2pi * i / F->Re()->GetUnitSteps();
		for (int j = 0; j < _N; j++)
		{
			_re += data[j] * cos((double)_2pi_f*(_t0+j/_UnitSteps));
			_im += -data[j] * sin((double)_2pi_f*(_t0 + j / _UnitSteps));
		}

		F->Re()->SetAt(i, _re/_UnitSteps);
		F->Im()->SetAt(i, _im / _UnitSteps);
	}

	return F;
}



//create n-th derivate function of this function 
Function * Function::CalculateDerivate(unsigned int n){
	
	if (!data)
		return NULL;

	Function* f = this;
	if (n>1)
		f = CalculateDerivate(n - 1);


	Function* df = new Function(f);

	for (long long i = 0; i < f->_N - 1; i++)
		df->SetAt(i, ((*f)[i + 1] - (*f)[i])*_UnitSteps);

	df->SetAt(_N - 1, 0);
	if (n>1) //delete temp n-1 order derivative
		delete f;

	return df;
}

void Function::ComputePartialDerivative(int k,Function* df, Function* f, Function* f1,double unitSteps)
{
	if (k == 0)
		unitSteps = f->_UnitSteps;

	for (long long i = 0; i < f->_N; i++) {
		if (f->fdata){
			if (k == 0){
				if (i < f->_N - 1) {
					ComputePartialDerivative(-1, df->fdata[i], f->fdata[i], f->fdata[i + 1], unitSteps);
				}
				else {
					ComputePartialDerivative(-1, df->fdata[i], f->fdata[i], NULL, unitSteps);
				}
			}
			else{
				ComputePartialDerivative(k-1, df->fdata[i], f->fdata[i], 
					f1 ? f1->fdata[i] : NULL, unitSteps);
			}
		}
		else{
			if (k == 0 && i < f->_N - 1) {
				df->data[i] = (f->data[i + 1] - f->data[i])*unitSteps;
			}
			else if (f1){
				df->data[i] = (f1->data[i] - f->data[i])*unitSteps;
			}
			else
				df->data[i] = 0;
		}
	}
}

Function * Function::CalculatePartialDerivate(unsigned int k, unsigned int n)
{
	Function* f = this;
	if (n>1)
		f = CalculatePartialDerivate(k,n - 1);

	Function* df = new Function(f);
	ComputePartialDerivative(k,df, f);

	if (n>1) //delete temp n-1 order derivative
		delete f;

	return df;
	
}

//volume integration of multi args F({x0,..xn-1}) on {[a0,b0],...,[an-1,bn-1]}
double Function::Integrate(double* a, double* b,int N)
{
	//use monte-carlo method
	std::random_device rd;
	std::mt19937 rgen(rd());
	std::uniform_real_distribution<> rdistr;

	double V = 1;
	for (int i = 0; i < _level; i++) {
		V *= (b[i] - a[i]);
	}
	double sum = 0;
	double* X = new double[_level];
	for (int i = 0; i < N; i++) {
		for (int i = 0; i < _level; i++) {
			X[i] = a[i] + (b[i] - a[i]) * rdistr(rgen);
		}
		sum += (*this)[X];
	}
	delete X;
	return V/N*sum;
}

//integrate on the [t1,t2] interval
double Function::Integrate(double t1, double t2)
{
	if (!data)
		return 0;

	if (t1 > t2)
		return Integrate(t2, t1);

	if (t2 > _T || t1 < _t0 || t1 == t2)
		return 0;

	double sum = 0;
	double u2 = 2 * _UnitSteps;
	long long sn = (long long)((t1 - _t0)*_UnitSteps);
	long long en = (long long)(_N - (_T - t2)*_UnitSteps);
	for (long long i = sn; i<en - 1; i++)
		sum += (data[i] + data[i + 1]) / u2;

	return sum;
}

//normalize the function -- divide by sqrt((f*f)->Integral())
void Function::Normalize() {

	//calculate f*f integral
	double sum = 0;
	double u2 = 2 * _UnitSteps;
	long long computeSteps = GetSize() - 1;
	for (long long i = 0; i < computeSteps; i++)
		sum += (pow(data[i],2) + pow(data[i + 1],2)) / u2;

	double norm = sqrt(sum); 

	for (long long i = 0; i < computeSteps; i++)
		data[i] /= norm;
}

void Function::Reset(double val) {
	if (!data) return;

	if (val == 0)
		memset(data, 0x00, _N*sizeof(double));
	else {
		for (long long i = 0; i < _N; i++)
			data[i] = val;
	}



}
void Function::ConvergeRight(double TMax) {
	if (!data) return;
	int n = (int)((TMax - _t0)*_UnitSteps);
	if (n > 0 && n < _N - 1) {
		double val = data[n];
		for (int i = n + 1; i < _N; i++){
			data[i] = val*exp(double(n-i)/(_N-i));
		}
	}
}
//use to set T to TMax
long long Function::CutRight(double TMax) {
	if (!data) return 0;
	long long oldN = _N;
	int n = (int)((TMax - _t0)*_UnitSteps);
	if (n > 0 && n < _N - 1) {
#pragma omp critical
		{maxN += (_N - n);}
		
		_N = n;
		_T = TMax;
	}
	return oldN;
}
//inverse of CutRight _N to old one returned by CutRight
//use the value returned by CutRight as parameter oldN
void Function::UndoCutRight(long long oldN) {
	if (!data) return;
	if (oldN > _N) {
#pragma omp critical
		{ maxN -= (oldN - _N); }
		_N = oldN;
		_T = _t0 + (_N - 1)/_UnitSteps;
	}
}

//create a function I(t) where I(t) is integral over[t0,t] of this function
Function * Function::CalculateIntegral(){
	if (!data)
		return NULL;

	Function *f = new Function(this);
	long long computeSteps = GetSize() - 1;
	double u2 = 2 * _UnitSteps;
	double t = _t0, sum = 0;
	f->SetAt(0, 0);
	for (long long i = 0; i < computeSteps; i++)
	{
		sum += (data[i] + data[i + 1])/u2;
		f->SetAt(i + 1, sum);
	}

	return f;
}
Function* Function::operator*(const Function* p) {
	if (p->_T != _T || p->_t0 != _t0 || p->_UnitSteps != _UnitSteps)
		return NULL;

	if (!data) {
		return NULL;//TODO implement operator for multi args functions
	}

	Function* f = new Function(this);
	for (long long i = 0; i < _N; i++)
		f->SetAt(i, data[i] * (p->data[i]));

	return f;

}

Function* Function::operator*(double x) {
	if (!data) {
		return NULL;//TODO implement operator for multi args functions
	}

	Function* f = new Function(this);
	for (long long i = 0; i < _N; i++)
		f->SetAt(i, data[i] * x);

	return f;
}

Function* Function::operator/(const Function* p) {
	if (p->_T != _T || p->_t0 != _t0 || p->_UnitSteps != _UnitSteps)
		return NULL;

	if (!data) {
		return NULL;//TODO implement operator for multi args functions
	}

	Function* f = new Function(this);
	for (long long i = 0; i < _N; i++)
		f->SetAt(i, p->data[i] == 0 ? HUGE_VAL : data[i] / (p->data[i]));

	return f;

}

Function* Function::operator/(double x) {
	if (!data) {
		return NULL;//TODO implement operator for multi args functions
	}

	Function* f = new Function(this);
	for (long long i = 0; i < _N; i++)
		f->SetAt(i, data[i] / x);

	return f;
}

Function* Function::operator+(const Function* p) {
	if (p->_T != _T || p->_t0 != _t0 || p->_UnitSteps != _UnitSteps)
		return NULL;

	if (!data) {
		return NULL;//TODO implement operator for multi args functions
	}


	Function* f = new Function(this);

	for (long long i = 0; i < _N; i++)
		f->SetAt(i, data[i] + (p->data[i]));

	return f;

}

Function* Function::operator+(double x) {
	if (!data) {
		return NULL;//TODO implement operator for multi args functions
	}

	Function* f = new Function(this);
	for (long long i = 0; i < _N; i++)
		f->SetAt(i, data[i] + x);

	return f;
}

Function* Function::operator-(const Function* p) {
	if (p->_T != _T || p->_t0 != _t0 || p->_UnitSteps != _UnitSteps)
		return NULL;

	if (!data) {
		return NULL;//TODO implement operator for multi args functions
	}

	Function* f = new Function(this);
	for (long long i = 0; i < _N; i++)
		f->SetAt(i, data[i] - (p->data[i]));

	return f;

}

Function* Function::operator-(double x) {
	if (!data) {
		return NULL;//TODO implement operator for multi args functions
	}

	Function* f = new Function(this);
	for (long long i = 0; i < _N; i++)
		f->SetAt(i, data[i] + x);

	return f;
}

Function* Function::operator-() {
	if (!data) {
		return NULL;//TODO implement operator for multi args functions
	}

	Function* f = new Function(this);
	for (long long i = 0; i < _N; i++)
		f->SetAt(i, -data[i]);

	return f;
}
Function* Function::SubstractFrom(double x){
	if (!data) {
		return NULL;//TODO implement operator for multi args functions
	}

	Function* f = new Function(this);
	for (long long i = 0; i < _N; i++)
		f->SetAt(i, x - data[i] );

	return f;
}
//Analyze min,max,nodes -- print result on stdout

int Function::maxNOutputLines = 50;
void Function::Analyze(AnalyzesResult& result, FILE* file, bool bReset){

	if (!data) return;

	double _prevdoubles[2] = { data[0], data[1] };

	if (bReset)
	{
		result._MinVal = data[0];
		result._min = _t0;
		result._max = _t0;
		result._MaxVal = data[0];
		result.nodes = 0;
		result.localMins = 0;
		result.localMaxs = 0;

	}

	result._ConvergenceTime = _t0;
	result._ConvergenceDelta = 0;
	double _ConvergenceStart = _t0, _ConvergenceEnd = _t0, _Convergencedouble = fabs(data[1]);

	double lastMaxVal = -HUGE_VAL;
	double lastMaxTime, lastMinTime, lastNodeTime, nodePeriod = 0, nodeTime;
	double lastMaxPeriod, lastMinPeriod;
	int maxTimes = -1, minTimes = -1, nodeTimes = -1;
	double lastMinVal = HUGE_VAL;
	int OutputLines = 0;
	for (int i = 2; i < _N; i++)
	{
		if (_prevdoubles[1] == data[i])
			continue;

		double pt = _t0 + (i - 1)/_UnitSteps; // previous time
		double t = _t0 + i / _UnitSteps; // current time

		if (data[i] == 0 || (data[i] * _prevdoubles[1] < 0 && data[i] * _prevdoubles[1] > -100 / _UnitSteps /*avoid divergence detection as nodes, for example tangent*/)) // node detected
		{
			if (result.nodesTimeSize > result.nodes)
				result.nodesTime[result.nodes] = t;

			result.nodes++;
			_Convergencedouble = 0;// reset convergence calc

			if (nodePeriod == 0 && nodeTimes == 0)
			{
				nodePeriod = t - lastNodeTime;
				nodeTimes = 1;
				lastNodeTime = t;
			}
			else if (nodeTimes > 0 && fabs(nodePeriod - t + lastNodeTime) <= 2 / _UnitSteps) //period detected
			{
				lastNodeTime = t;
				nodeTimes++;
			}
			else
			{
				if (nodeTimes > 1)
				{
					if (file && ++OutputLines < maxNOutputLines) {
						fprintf(file, "node @ t\t" double_FORMAT, nodeTime);
						fprintf(file, "\twith period \t" double_FORMAT " for %d times\n", nodePeriod, nodeTimes + 1);
					}
				}
				else if (nodeTimes == 1)
				{
					if (file && ++OutputLines < maxNOutputLines){
						fprintf(file, "node @ t\t" double_FORMAT "\n", nodeTime);
						fprintf(file, "node @ t\t" double_FORMAT "\n", lastNodeTime);
					}
				}
				nodeTimes = 0;
				nodePeriod = 0;
				lastNodeTime = t;
				nodeTime = t;
			}

			//if(file) fprintf(file,"\t\t[DBG]node @ t\t" double_FORMAT " count:%d\n",t,nodeTimes);
		}
		else
		{
			if (fabs(data[i]) < _Convergencedouble)
			{
				_Convergencedouble = fabs(data[i]);
				_ConvergenceEnd = t;
				if (_ConvergenceEnd - _ConvergenceStart > result._ConvergenceDelta)
				{
					result._ConvergenceDelta = _ConvergenceEnd - _ConvergenceStart;
					result._ConvergenceTime = t;
				}
			}
			else // reset convergence calc
			{
				_Convergencedouble = fabs(data[i]);
				_ConvergenceStart = _ConvergenceEnd = t;
			}
		}

		if (data[i] < _prevdoubles[1] && _prevdoubles[1] > _prevdoubles[0])
		{
			result.localMaxs++;
			if (fabs(lastMaxVal - _prevdoubles[1]) > 2 / _UnitSteps)
			{
				if (file && maxTimes >= 0 && ++OutputLines < maxNOutputLines){
					fprintf(file, "local max @ t\t" double_FORMAT "\ty[t]:\t" double_FORMAT, lastMaxTime, lastMaxVal);
					if (maxTimes > 0)
						fprintf(file, "\twith period \t" double_FORMAT " for %d times\n", lastMaxPeriod, maxTimes);
					else
						fprintf(file, "\n");
				}

				maxTimes = 0;
				lastMaxVal = _prevdoubles[1];
				lastMaxTime = pt;
				lastMaxPeriod = 0;
			}
			else //period detected
			{
				if (lastMaxPeriod == 0)
					lastMaxPeriod = pt - lastMaxTime;
				maxTimes++;
			}
			//if(file) fprintf(file,"\t\t[DBG]max @ t\t" double_FORMAT " count:%d\n",pt,maxTimes);
		}

		if (data[i] > _prevdoubles[1] && _prevdoubles[1] < _prevdoubles[0])
		{
			result.localMins++;
			if (fabs(lastMinVal - _prevdoubles[1]) > 2 / _UnitSteps)
			{
				if (file && minTimes >= 0 && ++OutputLines < maxNOutputLines){
					fprintf(file, "local min @ t\t" double_FORMAT "\ty[t]:\t" double_FORMAT, lastMinTime, lastMinVal);
					if (minTimes > 0)
						fprintf(file, "\twith period\t" double_FORMAT " for %d times\n", lastMinPeriod, minTimes);
					else
						fprintf(file, "\n");
				}

				lastMinVal = _prevdoubles[1];
				minTimes = 0;
				lastMinTime = pt;
				lastMinPeriod = 0;
			}
			else
			{
				if (lastMinPeriod == 0)
					lastMinPeriod = pt - lastMinTime;
				minTimes++;
			}
			//if(file) fprintf(file,"\t\t[DBG]min @ t\t" double_FORMAT " count:%d\n",pt,minTimes);
		}

		if (result._MaxVal < data[i])
		{
			result._max = t;
			result._MaxVal = data[i];
		}
		if (result._MinVal > data[i]) {
			result._min = t;
			result._MinVal = data[i];
		}

		_prevdoubles[0] = _prevdoubles[1];
		_prevdoubles[1] = data[i];
	}

	if (file)
	{
		if (OutputLines < maxNOutputLines) {
			if (nodeTimes > 1) {
				fprintf(file, "node @ t\t" double_FORMAT, nodeTime);
				fprintf(file, "\twith period \t" double_FORMAT " for %d times\n", nodePeriod, nodeTimes);
			}
			else if (nodeTimes == 1){
				fprintf(file, "node @ t\t" double_FORMAT "\n", nodeTime);
				fprintf(file, "node @ t\t" double_FORMAT "\n", lastNodeTime);

			}
			else if (nodeTimes == 0){
				fprintf(file, "node @ t\t" double_FORMAT "\n", nodeTime);
			}

			if (maxTimes >= 0){
				fprintf(file, "local max @ t\t" double_FORMAT "\ty[t]:\t" double_FORMAT, lastMaxTime, lastMaxVal);
				if (maxTimes > 0)
					fprintf(file, "\twith period \t" double_FORMAT " for %d times\n", lastMaxPeriod, maxTimes);
				else
					fprintf(file, "\n");
			}

			if (minTimes >= 0){
				fprintf(file, "local min @ t\t" double_FORMAT "\ty[t]:\t" double_FORMAT, lastMinTime, lastMinVal);
				if (minTimes > 0)
					fprintf(file, "\twith period\t" double_FORMAT " for %d times\n", lastMinPeriod, minTimes);
				else
					fprintf(file, "\n");
			}
		}
		else {
			fprintf(file, "........\nWARNING: Too many lines. Output was cut.\n");
		}

		fprintf(file, "-----\nglobal min @ t\t" double_FORMAT "\ty[t]:\t" double_FORMAT "\n", result._min, result._MinVal);
		fprintf(file, "global max @ t\t" double_FORMAT "\ty[t]:\t" double_FORMAT "\n", result._max, result._MaxVal);
		fprintf(file, "convergence @ t\t" double_FORMAT "\tdelta:\t" double_FORMAT "\n", result._ConvergenceTime, result._ConvergenceDelta);
		fprintf(file, "%d nodes \n%d local min(s)\n%d local max(s)\n-----\n", result.nodes, result.localMins, result.localMaxs);


	}
}

//Dump all data one per row to console
void Function::Dump(FILE* file){
	for (int i = 0; i < _N; i++)
		fprintf(file, double_FORMAT "\n", data[i]);
}


//1st parameter of op is current double of t
//2nd parameter of op is this function double, use only doubles of function before or including step
//3rd paramter of op is current calculations step
Function* Function::ApplyOperator(const Operator& O){
	if (!data) return NULL;

	Function* f = new Function(this);
	long long computeSteps = GetSize();


	for (long long j = 0; j < computeSteps; j++)
		f->SetAt(j, O.op(_t0 + j / _UnitSteps, *this, j));

	return f;
}

ComplexFunction* Function::ApplyOperator(const ComplexOperator<Function>& O){
	
	if (!data) return NULL;

	ComplexFunction* f = new ComplexFunction(_T, 0, 0, _t0, GetPrecision());
	long long computeSteps = GetSize();
	

	for (long long j = 0; j < computeSteps; j++)
	{
		double t = _t0 + j / _UnitSteps;
		f->Re()->SetAt(j, O.opRe(t, *this, j));
		f->Im()->SetAt(j, O.opIm(t, *this, j));
	}

	return f;
}

Function::operator ComplexFunction*() {
	if (!data) return NULL;

	Function* Re = new Function(this, true);
	Function* Im = new Function(this);
	Im->ResetData();

	return new ComplexFunction(Re, Im);
}

Function* Function::ApplyOperatorNL(OperatorDef op){
	if (!data) return NULL;

	Function* f = new Function(this);
	long long computeSteps = GetSize();
	

	for (long long j = 0; j < computeSteps; j++)
		f->SetAt(j, op(_t0 + j / _UnitSteps, *this, j));

	return f;
}

ComplexFunction* Function::ApplyComplexOperatorNL(OperatorDef op_Re, OperatorDef op_Im){
	if (!data) return NULL;

	ComplexFunction* f = new ComplexFunction(this);
	long long computeSteps = GetSize();
	

	for (long long j = 0; j < computeSteps; j++)
	{
		double t = _t0 + j / _UnitSteps;
		f->Re()->SetAt(j, op_Re(t, *this, j));
		f->Im()->SetAt(j, op_Im(t, *this, j));
	}

	return f;
}


ComplexFunction::ComplexFunction(double T, double Re0, double Im0, double t0, int precision){
	_Re = new Function(T, Re0, t0, precision);
	_Im = new Function(T, Im0, t0, precision);
}

ComplexFunction::ComplexFunction(Function* f){
	_Re = new Function(f);
	_Im = new Function(f);
}

ComplexFunction::ComplexFunction(Function* real, Function* imag)
{
	if (!real) {
		_Re = new Function(imag);
		_Re->Reset();
	}
	else
		_Re = real;
	
	if (!imag) {
		_Im = new Function(real);
		_Im->Reset();
	}
	else
		_Im = imag;
}

ComplexFunction::ComplexFunction(double T, std::function<double(double)> genRe, std::function<double(double)> genIm,
	double t0, int precision){
	_Re = new Function(T, genRe, t0, precision);
	_Im = new Function(T, genIm, t0, precision);
}

ComplexFunction::ComplexFunction(const ComplexFunction *f, bool bCopyData)
{
	_Re = new Function(f->_Re,bCopyData);
	_Im = new Function(f->_Im, bCopyData);
}

ComplexFunction::~ComplexFunction(){
	delete _Re;
	delete _Im;
}



Function* ComplexFunction::InverseFourierTransform(double T, double t0, int precision){
	Function* invF = new Function(T, 0.0, t0, precision);
	long long computesteps = invF->GetSize();
	const double _2pi = 2 * π;
	long long N = _Re->GetSize();

#pragma omp parallel for
	for (long long i = 0; i < computesteps; i++)
	{
		double sum = 0;
		double _2pi_t = _2pi * (t0 + i / invF->GetUnitSteps());
		for (long long j = 0; j < N; j++){
			double _2pi_t_f = _2pi_t * j / _Re->GetUnitSteps();
			sum += (*_Re)[j] * cos(_2pi_t_f) - (*_Im)[j] * sin(_2pi_t_f);
		}
		invF->SetAt(i, sum / _Re->GetUnitSteps() * 2);
	}

	return invF;
}

Function* ComplexFunction::Module(){
	Function* m = new Function(_Re);
	long long computesteps = _Re->GetSize();

	for (long long i = 0; i<computesteps; i++)
		m->SetAt(i, sqrt(pow((*_Re)[i], 2) + pow((*_Im)[i], 2)));
	return m;
}

Function* ComplexFunction::Arg(){
	Function* arg = new Function(_Re);
	long long computesteps = _Re->GetSize();

	for (long long i = 0; i<computesteps; i++)
		arg->SetAt(i, atan2(_Im->data[i],_Re->data[i]));
	return arg;
}

//normalize the function -- divide by sqrt((f*f)->Integral())
void ComplexFunction::Normalize() {
	double sum = pow(_Re->data[0], 2) + pow(_Im->data[0], 2) / 2;
	long long C = _Re->GetSize() - 1;
	for (long long i = 1; i < C; i++)
		sum += (pow(_Re->data[i], 2) + pow(_Im->data[i], 2));

	sum += pow(_Re->data[C], 2) + pow(_Im->data[C], 2) / 2;
	double norm = sqrt(sum / _Re->_UnitSteps);

	for (long long i = 0; i <= C; i++){
		_Re->data[i] /= norm;
		_Im->data[i] /= norm;
	}

}

ComplexFunction* ComplexFunction::operator*(const ComplexFunction* p) {
	if (p->_Re->_T != _Re->_T || p->_Re->_t0 != _Re->_t0 || p->_Re->_UnitSteps != _Re->_UnitSteps)
		return NULL;

	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] * (p->_Re->data[i]) - _Im->data[i] * (p->_Im->data[i]));
		f->_Im->SetAt(i, _Re->data[i] * (p->_Im->data[i]) + _Im->data[i] * (p->_Re->data[i]));
	}

	return f;

}

ComplexFunction* ComplexFunction::operator*(const Function* p) {
	if (p->_T != _Re->_T || p->_t0 != _Re->_t0 || p->_UnitSteps != _Re->_UnitSteps)
		return NULL;

	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] * (p->data[i]));
		f->_Im->SetAt(i, _Im->data[i] * (p->data[i]));
	}

	return f;

}

ComplexFunction* ComplexFunction::operator*(double x) {
	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] * x);
		f->_Im->SetAt(i, _Im->data[i] * x);
	}

	return f;
}
ComplexFunction* ComplexFunction::operator*(const Complex& x){
	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] * x.real() - _Im->data[i] * x.imag());
		f->_Im->SetAt(i, _Re->data[i] * x.imag() + _Im->data[i] * x.real());
	}

	return f;
}

//this z1/z2 = z1*(~z2)/ (|z2|^2)
ComplexFunction* ComplexFunction::operator/(const ComplexFunction* p) {
	if (p->_Re->_T != _Re->_T || p->_Re->_t0 != _Re->_t0 || p->_Re->_UnitSteps != _Re->_UnitSteps)
		return NULL;

	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		double denom = pow(p->_Re->data[i], 2) + pow(p->_Im->data[i], 2); //|z2|^2
		if (denom == 0) {
			f->_Re->SetAt(i, HUGE_VAL);
			f->_Im->SetAt(i, HUGE_VAL);
		}
		else {
			f->_Re->SetAt(i, (_Re->data[i] * (p->_Re->data[i]) + _Im->data[i] * (p->_Im->data[i])) / denom);
			f->_Im->SetAt(i, (_Im->data[i] * (p->_Re->data[i]) - _Re->data[i] * (p->_Im->data[i])) / denom);
		}
	}

	return f;

}

ComplexFunction* ComplexFunction::operator/(const Function* p) {
	if (p->_T != _Re->_T || p->_t0 != _Re->_t0 || p->_UnitSteps != _Re->_UnitSteps)
		return NULL;

	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		if (p->data[i] == 0){
			f->_Re->SetAt(i, HUGE_VAL);
			f->_Im->SetAt(i, HUGE_VAL);
		}
		else {
			f->_Re->SetAt(i, _Re->data[i] / (p->data[i]));
			f->_Im->SetAt(i, _Im->data[i] / (p->data[i]));
		}
	}

	return f;

}

ComplexFunction* ComplexFunction::operator/(double x) {
	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] / x);
		f->_Im->SetAt(i, _Im->data[i] / x);
	}

	return f;
}

//this z1/z2 = z1*(~z2)/ (|z2|^2)
ComplexFunction* ComplexFunction::operator/(const Complex& x){
	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		double denom = pow(x.real(), 2) + pow(x.imag(), 2); //|z2|^2
		if (denom == 0) {
			f->_Re->SetAt(i, HUGE_VAL);
			f->_Im->SetAt(i, HUGE_VAL);
		}
		else {
			f->_Re->SetAt(i, (_Re->data[i] * x.real() + _Im->data[i] * x.imag()) / denom);
			f->_Im->SetAt(i, (_Im->data[i] * x.real() - _Re->data[i] * x.imag()) / denom);
		}
	}

	return f;
}
//scalar product (inner product)
Complex ComplexFunction::operator | (const ComplexFunction* p){
	if (p->_Re->_T != _Re->_T || p->_Re->_t0 != _Re->_t0 || p->_Re->_UnitSteps != _Re->_UnitSteps)
		return NULL;

	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	long long i = 0;
	Complex res( (_Re->data[i] * (p->_Re->data[i]) + _Im->data[i] * (p->_Im->data[i])) / 2,
				 (_Re->data[i] * (p->_Im->data[i]) - _Im->data[i] * (p->_Re->data[i])) / 2) ;
	for (i = 1; i < _Re->_N - 1; i++) {
		res._Val[0] += _Re->data[i] * (p->_Re->data[i]) + _Im->data[i] * (p->_Im->data[i]);
		res._Val[1] += _Re->data[i] * (p->_Im->data[i]) - _Im->data[i] * (p->_Re->data[i]);
	}

	res._Val[0] += (_Re->data[i] * (p->_Re->data[i]) + _Im->data[i] * (p->_Im->data[i])) / 2;
	res._Val[1] += (_Re->data[i] * (p->_Im->data[i]) - _Im->data[i] * (p->_Re->data[i])) / 2;

	res._Val[0] /= _Re->_UnitSteps;
	res._Val[1] /= _Re->_UnitSteps;

	return res;
}

ComplexFunction* ComplexFunction::operator+(const ComplexFunction* p) {
	if (p->_Re->_T != _Re->_T || p->_Re->_t0 != _Re->_t0 || p->_Re->_UnitSteps != _Re->_UnitSteps)
		return NULL;

	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] + (p->_Re->data[i]));
		f->_Im->SetAt(i, _Im->data[i] + (p->_Im->data[i]));
	}

	return f;

}

ComplexFunction* ComplexFunction::operator+(const Function* p) {
	if (p->_T != _Re->_T || p->_t0 != _Re->_t0 || p->_UnitSteps != _Re->_UnitSteps)
		return NULL;

	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] + (p->data[i]));
		f->_Im->SetAt(i, _Im->data[i]);
	}

	return f;

}

ComplexFunction* ComplexFunction::operator+(double x) {
	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (int i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] + x);
		f->_Im->SetAt(i, _Im->data[i]);
	}

	return f;
}

ComplexFunction* ComplexFunction::operator+(const Complex& x){
	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (int i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] + x.real());
		f->_Im->SetAt(i, _Im->data[i] + x.imag());
	}

	return f;
}

ComplexFunction* ComplexFunction::operator-(const ComplexFunction* p) {
	if (p->_Re->_T != _Re->_T || p->_Re->_t0 != _Re->_t0 || p->_Re->_UnitSteps != _Re->_UnitSteps)
		return NULL;

	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] - (p->_Re->data[i]));
		f->_Im->SetAt(i, _Im->data[i] - (p->_Im->data[i]));
	}

	return f;

}

ComplexFunction* ComplexFunction::operator-(const Function* p) {
	if (p->_T != _Re->_T || p->_t0 != _Re->_t0 || p->_UnitSteps != _Re->_UnitSteps)
		return NULL;

	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (long long i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] - (p->data[i]));
		f->_Im->SetAt(i, _Im->data[i]);
	}

	return f;

}

ComplexFunction* ComplexFunction::operator-(double x) {
	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (int i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] + x);
		f->_Im->SetAt(i, _Im->data[i]);
	}

	return f;
}

ComplexFunction* ComplexFunction::operator-(const Complex& x){
	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (int i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i] - x.real());
		f->_Im->SetAt(i, _Im->data[i] - x.imag());
	}

	return f;
}



ComplexFunction* ComplexFunction::ScalarMinusFunction(double x, ComplexFunction& f1){
	if (!f1._Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(&f1);

	for (int i = 0; i < f1._Re->_N; i++) {
		f->_Re->SetAt(i, x - f1._Re->data[i]);
		f->_Im->SetAt(i, f1._Im->data[i]);
	}

	return f;
}

ComplexFunction* ComplexFunction::ScalarMinusFunction(const Complex& x, ComplexFunction& f1){
	if (!f1._Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(&f1);

	for (int i = 0; i < f1._Re->_N; i++) {
		f->_Re->SetAt(i, x.real() - f1._Re->data[i]);
		f->_Im->SetAt(i, x.imag() - f1._Im->data[i]);
	}

	return f;
}

ComplexFunction* ComplexFunction::ScalarMinusFunction(const Function& x, ComplexFunction& f1){
	if (!f1._Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(&f1);

	for (int i = 0; i < f1._Re->_N; i++) {
		f->_Re->SetAt(i, x.data[i] - f1._Re->data[i]);
		f->_Im->SetAt(i, f1._Im->data[i]);
	}

	return f;
}
ComplexFunction* ComplexFunction::operator-() {
	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (int i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, -_Re->data[i]);
		f->_Im->SetAt(i, -_Im->data[i]);
	}

	return f;
}

ComplexFunction* ComplexFunction::operator~(){
	if (!_Re->data) {
		return NULL;//TODO implement operator for multi args functions
	}

	ComplexFunction* f = new ComplexFunction(this);

	for (int i = 0; i < _Re->_N; i++) {
		f->_Re->SetAt(i, _Re->data[i]);
		f->_Im->SetAt(i, -_Im->data[i]);
	}

	return f;
}

//1st parameter of op is current double of t
//2nd parameter of op is this function double, use only doubles of function before or including step
//3rd paramter of op is current calculations step
//use NL version to reduce process time
ComplexFunction* ComplexFunction::ApplyOperator(const ComplexOperator<>& O)
{
	ComplexFunction* f = new ComplexFunction(_Re->GetT(), 0, 0, _Re->GetT0(), _Re->GetPrecision());
	long long computeSteps = _Re->GetSize();

	for (long long j = 0; j < computeSteps; j++)
	{
		double t = _Re->GetT0() + j / _Re->GetUnitSteps();
		f->Re()->SetAt(j, O.opRe(t, *this, j));
		f->Im()->SetAt(j, O.opIm(t, *this, j));
	}

	return f;
}

ComplexFunction* ComplexFunction::ApplyComplexOperatorNL(OperatorDef op_Re, OperatorDef op_Im){
	ComplexFunction* f = new ComplexFunction(_Re->GetT(), 0, 0, _Re->GetT0(), _Re->GetPrecision());
	long long computeSteps = _Re->GetSize();


	for (long long j = 0; j < computeSteps; j++)
	{
		double t = _Re->GetT0() + j / _Re->GetUnitSteps();
		f->Re()->SetAt(j, op_Re(t, *this, j));
		f->Im()->SetAt(j, op_Im(t, *this, j));
	}

	return f;
}

FunctionSystem::FunctionSystem(int n, double T, double* initialdoubles, double t0, int precision){
	_n = n;
	_Y = (Function **)malloc(n*sizeof(void*));
	for (int i = 0; i < _n; i++) {
		_Y[i] = new Function(T, initialdoubles[i], t0, precision);
	}
}

FunctionSystem::FunctionSystem(int n, int m, double* XMax, double *XMin, std::function<double(double* X, int k)> gen,int* precision) {
	_n = n;
	_Y = (Function **)malloc(n*sizeof(void*));
	for (int i = 0; i < _n; i++) {
		_Y[i] = new Function(m, XMax, XMin, 
			[i, gen](double* X) { 
				return gen(X, i); 
			},precision);
	}
}

FunctionSystem::~FunctionSystem(void){
	for (int i = 0; i < _n; i++)
		if (_Y[i]) delete _Y[i];
	free(_Y);
}

void FunctionSystem::Init(ODEAnalyzesResult*res, double t0, double initialVal){
	res->_MaxVal = res->_MinVal = initialVal;
	res->_min = t0;
	res->_max = t0;
	res->_MaxVal = initialVal;
	res->nodes = 0;
	res->_ConvergenceTime = t0;
	res->_ConvergenceDelta = 0;
	res->_ConvergenceStart = res->_ConvergenceEnd = t0;
	res->_Convergencedouble = fabs(initialVal);
	res->_prevdouble = initialVal;
}

void FunctionSystem::Compute(ODEAnalyzesResult*res, double t, double val){

	if (res->_prevdouble == val)
		return;

	if (val == 0 || val*res->_prevdouble < 0) // node detected
	{
		res->nodes++;
		res->_Convergencedouble = 0;// reset convergence calc
	}
	else
	{
		if (fabs(val) < res->_Convergencedouble)
		{
			res->_Convergencedouble = fabs(val);
			res->_ConvergenceEnd = t;
			if (res->_ConvergenceEnd - res->_ConvergenceStart > res->_ConvergenceDelta)
			{
				res->_ConvergenceDelta = res->_ConvergenceEnd - res->_ConvergenceStart;
				res->_ConvergenceTime = t;
			}
		}
		else // reset convergence calc
		{
			res->_Convergencedouble = fabs(val);
			res->_ConvergenceStart = res->_ConvergenceEnd = t;
		}
	}

	if (res->_MaxVal < val)
	{
		res->_max = t;
		res->_MaxVal = val;
	}
	if (res->_MinVal > val) {
		res->_min = t;
		res->_MinVal = val;
	}

	res->_prevdouble = val;
}

Function* FunctionSystem::CalculateFunction(std::function<double(double, FunctionSystem& F)> gen){
	Function* f = new Function(_Y[0]);
	f->ResetLastIndex();
	for (int i = 0; i<_n; i++)
		_Y[i]->ResetLastIndex();
	long long computeSteps = _Y[0]->GetSize();

	for (int j = 0; j < computeSteps; j++)
	{
		f->SetAt(j, gen(_Y[0]->GetT0() + j / _Y[0]->GetUnitSteps(), *this));
		for (int i = 0; i<_n; i++)
			_Y[i]->PopLastValue();
	}

	return f;
}

Function * FunctionSystem::CalculateFunctionNL(FunctionDef gen, void* userdata){
	Function* f = new Function(_Y[0]);
	f->ResetLastIndex();
	for (int i = 0; i<_n; i++)
		_Y[i]->ResetLastIndex();
	long long computeSteps = _Y[0]->GetSize();

	for (int j = 0; j < computeSteps; j++)
	{
		f->SetAt(j, gen(_Y[0]->GetT0() + j / _Y[0]->GetUnitSteps(), *this, userdata));
		for (int i = 0; i<_n; i++)
			_Y[i]->PopLastValue();
	}

	return f;
}

/*Solve eqn system Fk'=f[k](t,F0,F2,..Fn-1) for k=0..n-1
- f[k](t,F0,F2,..Fn-1) is defined as a set of k=0..n-1 FunctionDef(s)
- use F(k) for F0,F2,...Fn-1 inside each f[k] function*/
bool FunctionSystem::SolveODENL(FunctionDef* f, void* userdata, Method method, unsigned int accuracy, ODEAnalyzesResult* res){
	double t0 = _Y[0]->GetT0();
	double T = _Y[0]->GetT();
	int setDataFreq = (int)Function::__Pow10(accuracy);
	double unitSteps = _Y[0]->GetUnitSteps() * setDataFreq;
	long long computeSteps = long long(unitSteps*(T - t0));
	int setDataStep = 0;

	if (res) { //init analyzes
		for (int i = 0; i < _n; i++) {
			Init(res + i, t0, (*_Y[i])[0.]);
		}
	}

	switch (method)
	{
	case RK4:
	{
		double* y = new double[_n];
		double* k1 = new double[_n];
		double* k2 = new double[_n];
		double* k3 = new double[_n];
		double* k4 = new double[_n];

		for (int j = 0; j < computeSteps - 1; j++)
		{
			double t = t0 + j / unitSteps; 
			for (int i = 0; i < _n; i++)
			{
				y[i] = *_Y[i];//save current doubles
				k1[i] = f[i](t, *this, userdata);//compute k1
			}

			//compute k2
			for (int i = 0; i < _n; i++)
				*_Y[i] = y[i] + k1[i] / 2 / unitSteps; //update current doubles for k2

			t = t0 + (2*j+1) / unitSteps / 2; //t=t+h/2
			for (int i = 0; i < _n; i++)
				k2[i] = f[i](t, *this, userdata);

			//compute k3
			for (int i = 0; i < _n; i++)
				*_Y[i] = y[i] + k2[i] / 2 / unitSteps;

			for (int i = 0; i < _n; i++)
				k3[i] = f[i](t, *this, userdata);

			//compute k4
			for (int i = 0; i < _n; i++)
				*_Y[i] = y[i] + k3[i] / unitSteps;

			t = t0 + (j + 1) / unitSteps;//t=t+h
			for (int i = 0; i < _n; i++)
				k4[i] = f[i](t, *this, userdata);

			//compute final double
			for (int i = 0; i < _n; i++) {
				*_Y[i] = y[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6 / unitSteps;
				if (res) Compute(res + i, t, *_Y[i]);
			}

			if (++setDataStep == setDataFreq)
			{
				setDataStep = 0;
				for (int i = 0; i < _n; i++)
					_Y[i]->SaveCurrentValue();
			}

		}
		delete[] y;
		delete[] k1;
		delete[] k2;
		delete[] k3;
		delete[] k4;
		for (int i = 0; i < _n; i++)
			_Y[i]->SaveCurrentValue();
	}
		break;
	}
	return true;
}

bool FunctionSystem::SolveODE(std::function<double(double t, FunctionSystem& F)>* f, Method method, unsigned int accuracy, ODEAnalyzesResult* res){
	double t0 = _Y[0]->GetT0();
	double T = _Y[0]->GetT();
	int setDataFreq = (int)Function::__Pow10(accuracy);
	double unitSteps = _Y[0]->GetUnitSteps() * setDataFreq;
	long long computeSteps = long long(unitSteps*(T - t0));
	int setDataStep = 0;

	if (res) { //init analyzes
		for (int i = 0; i < _n; i++) {
			Init(res + i, t0, (*_Y[i])[0.]);
		}
	}

	switch (method)
	{
	case RK4:
	{
		double* y = new double[_n];
		double* k1 = new double[_n];
		double* k2 = new double[_n];
		double* k3 = new double[_n];
		double* k4 = new double[_n];

		for (int j = 0; j < computeSteps - 1; j++)
		{
			double t = t0 + j / unitSteps;
			for (int i = 0; i < _n; i++)
			{
				y[i] = *_Y[i];//save current doubles
				k1[i] = f[i](t, *this);//compute k1
			}

			//compute k2
			for (int i = 0; i < _n; i++)
				*_Y[i] = y[i] + k1[i] / 2 / unitSteps; //update current doubles for k2

			t = t0 + (2 * j + 1) / unitSteps / 2; //t=t+h/2
			for (int i = 0; i < _n; i++)
				k2[i] = f[i](t, *this);

			//compute k3
			for (int i = 0; i < _n; i++)
				*_Y[i] = y[i] + k2[i] / 2 / unitSteps;

			for (int i = 0; i < _n; i++)
				k3[i] = f[i](t, *this);

			//compute k4
			for (int i = 0; i < _n; i++)
				*_Y[i] = y[i] + k3[i] / unitSteps;

			t = t0 + (j + 1) / unitSteps;//t=t+h
			for (int i = 0; i < _n; i++)
				k4[i] = f[i](t, *this);

			//compute final double
			for (int i = 0; i < _n; i++) {
				*_Y[i] = y[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6 / unitSteps;
				if (res) Compute(res + i, t, *_Y[i]);
			}

			if (++setDataStep == setDataFreq)
			{
				setDataStep = 0;
				for (int i = 0; i < _n; i++)
					_Y[i]->SaveCurrentValue();
			}

		}
		delete[] y;
		delete[] k1;
		delete[] k2;
		delete[] k3;
		delete[] k4;
		for (int i = 0; i < _n; i++)
			_Y[i]->SaveCurrentValue();
	}
		break;
	}
	return true;
}

double Function::PDEData::getRandomVal(double& fval, double* fref, double unitSteps, double* fref2, double unitSteps2) // continous functions
{
	double val = fval;
	if (state == initial) {
		if (Function::IsNull(fval)) {
			boundary = true;
		}
		else if (fval == HUGE_VAL)
		{
			threshold += (double)1 / valSteps;
			unk1[unkGlobalIndex++] = &fval;
		}
	}
	else {
		int tid = GetOMPIndex();
		if (unk1[tid] == &fval) 
		{
			for (int i = tid - 1; i >= 0; i--)
				if (unk1[i] == &fval) { // already assigned previously
				return SetDuplicateOMPUnk(i);
				}

			double rnd;
			if (state == select) {
				if (fref){
					rnd = rDelta - 1 + (2 - rDelta) * rdistr(rgen);
					val = *fref + maxDelta * rnd / unitSteps;
					if (fref2) {
						rnd = rDelta - 1 + (2 - rDelta) * rdistr(rgen);
						val += *fref2 + maxDelta * rnd / unitSteps2;
						val /= 2;
					}
				}
				else
					val = valMin + (valMax - valMin) * rdistr(rgen);
			}
			else { //compute
				rnd = -1 + 2 * rdistr(rgen);//random between -1 and 1
				val = unk2[tid] + rnd / valSteps;
			}

			SetOMPUnk(val);
			
		}
	}
	
	return val;
}

double Function::PDEData::operator()(int i, int j) {
	long long* s = (state == initial) ? this->s : GetOMPStepVector();
	double f_s = getRandomVal(f[s]);
	s[i]++; //s= (s0,...si+1,...)
	double f_si = getRandomVal(f[s], &f_s, UnitSteps[i]);

	if (j == -1) { //first order partial derivative ∂f/∂xi
		s[i]--; //initial s
		return (f_si - f_s)*UnitSteps[i];
	}
	else if (i == j){ 
		//∂^2f/∂xi^2
		s[i]++;  //s= (s0,...si+2,....)
		double f_sii = getRandomVal(f[s], &f_si, UnitSteps[i]);
		s[i] -= 2; //initial s
		return (f_sii - 2 * f_si + f_s)*pow(UnitSteps[i], 2);
	
	}else{
		// ∂^2f/∂xi∂xj
		s[i]--; s[j]++;  //s= (s0,...si,...,sj+1,....)
		double f_sj = getRandomVal(f[s], &f_s, UnitSteps[j]);
		s[i]++;//s= (s0,...si+1,...,sj+1,....)
		double& ref_sij = f[s];
		double f_sij = getRandomVal(ref_sij, &f_si, UnitSteps[i], &f_sj, UnitSteps[j]);
		s[j]--; s[i]--; //initial s
		return (f_sij - f_si - f_sj + f_s)*UnitSteps[i] * UnitSteps[j];
	}
	
}


long long* Function::PDEData::GetOMPStepVector(){
	int thread_no = 0;
#ifdef _OPENMP
	thread_no = omp_get_thread_num();
#endif

	return omps + thread_no*n;
}

void Function::PDEData::InitOMPData() {
	if (!unkd) {//execute only first time
		int threads_no = 1;
#ifdef _OPENMP
		threads_no = omp_get_num_threads();
#endif

		unkd = new double[threads_no * unkMaxSize];
		unkIndex = new int[threads_no];
		omps = new long long[threads_no*n];
	}
}

double Function::PDEData::SetDuplicateOMPUnk(int i){
	int thread_no = 0;
#ifdef _OPENMP
	thread_no = omp_get_thread_num();
#endif
	double val = unkd[thread_no * unkMaxSize + i];
	unkd[thread_no * unkMaxSize + unkIndex[thread_no]] = val;
	unkIndex[thread_no]++;
	return val;
}

void Function::PDEData::SetOMPUnk(double val){
	int thread_no = 0;
#ifdef _OPENMP
	thread_no = omp_get_thread_num();
#endif

	unkd[thread_no * unkMaxSize + unkIndex[thread_no]] = val;
	unkIndex[thread_no]++;
}

void Function::PDEData::ResetOMPData(){
	int thread_no = 0;
#ifdef _OPENMP
	thread_no = omp_get_thread_num();
#endif
	unkIndex[thread_no] = 0;
	memcpy(omps + thread_no*n, s, n*sizeof(long long));
}

int Function::PDEData::GetOMPIndex(){
	int thread_no = 0;
#ifdef _OPENMP
	thread_no = omp_get_thread_num();
#endif
	return unkIndex[thread_no];
}

double* Function::PDEData::GetOMPUnk()
{
	int thread_no = 0;
#ifdef _OPENMP
	thread_no = omp_get_thread_num();
#endif

	return unkd + thread_no*unkMaxSize;

}
//pdeData.X is current vector {x0,x1,x2,...xn-1}
//pdeData.s is current discrete index vector {s0,s1,s2,...sn-1}
//pdeData.UnitSteps is the multiplicative inverse of delta vector {Δx0,Δx1,Δx2,...Δxn-1}
//l is current variable index , l start from 0
bool Function::SolvePDENL_MOLRecursive(int l, PDEData& pdeData) {
	long long steps = GetSize();
	pdeData.UnitSteps[l] = _UnitSteps;
	if (fdata) {
		
		for (long long s = 0; s < steps; s++) {
			pdeData.X[l] = _t0 + s / _UnitSteps;
			pdeData.s[l] = s;
			unsigned int tick;
			if (DBG_PDE_SOLVE && l == 0) {
				tick = os_ticks();
				pdeData.retriesCount = 0;
				pdeData.maxRetries = 0;
			}

			if (!fdata[s]->SolvePDENL_MOLRecursive(l + 1, pdeData))
				return false;
			
			if (DBG_PDE_SOLVE && l == 0) {
				printf("[%02d%%] %ums\tscore: %lg\tRetries:%02.2lf%% MaxRetries:%u\n", int(((s + 1) * 100) / steps), os_ticks() - tick, pdeData.score / (s + 1), pdeData.retriesCount,pdeData.maxRetries);
			}
		}

	}
	else {


		for (long long s = 0; s < steps; s++) {
			pdeData.X[l] = _t0 + s / _UnitSteps;
			pdeData.s[l] = s;
			pdeData.state = PDEData::initial;
			pdeData.unkGlobalIndex = 0;
			memset(pdeData.unk1, 0x00, pdeData.unkMaxSize*sizeof(double*));
			pdeData.threshold = 0;
			pdeData.state = PDEData::initial;
			pdeData.boundary = false;
			pdeData.eqn(pdeData.X, pdeData, pdeData.userdata);
			if (pdeData.unkGlobalIndex == 0)
				continue;//no unknowns - go next step
			if( pdeData.boundary)
			{
				//check at least current value is set
				double& fval = pdeData.f[pdeData.s];
				if (fval == HUGE_VAL) {// set same as previous value
					pdeData.s[l]--;
					fval = pdeData.f[pdeData.s];
				}
				continue; // bondary condition -- go next step
			}
			#pragma omp parallel
			{
			#pragma omp master
				{
					pdeData.InitOMPData();
				}
			}

			unsigned int tries = 0;
			double resMin;
			double initialDelta = pdeData.maxDelta;
			unsigned int selectSteps = pdeData.selectSteps;
			unsigned int computeSteps =  pdeData.computeSteps;

			do{
				resMin = HUGE_VAL;
				pdeData.state = PDEData::select;
				#pragma omp parallel //num_threads(1)
				{
					for (unsigned int c = 0; c < selectSteps; c++){
						pdeData.ResetOMPData();
						double res = abs(pdeData.eqn(pdeData.X, pdeData, pdeData.userdata));

						if (resMin > res) {
							#pragma omp critical
							{
								if (resMin > res) {
									resMin = res;
									//save current best values
									memcpy(pdeData.unk3, pdeData.GetOMPUnk(), pdeData.unkGlobalIndex * sizeof(double));

								}
							}
						}
					}
				}
					
				pdeData.state = PDEData::compute; // go into higher precision depth with selected candidate
				
				double intialValSteps = pdeData.valSteps;
				for (int depth = 0; depth < 10; depth++) {
					//set reference to current best values
					memcpy(pdeData.unk2, pdeData.unk3, pdeData.unkGlobalIndex * sizeof(double));

					#pragma omp parallel //num_threads(1)
					{
						for (unsigned int c = 0; c < computeSteps; c++){
							pdeData.ResetOMPData();
							double res = abs(pdeData.eqn(pdeData.X, pdeData, pdeData.userdata));

							if (resMin > res) {
								#pragma omp critical
								{
									if (resMin > res) {
										resMin = res;
										//save current best values
										memcpy(pdeData.unk3, pdeData.GetOMPUnk(), pdeData.unkGlobalIndex * sizeof(double));
									}
								}
							}
						}
					}
			
					pdeData.valSteps *= 4;
				}
				pdeData.valSteps = intialValSteps;
				//save current solution
				for (int p = 0; p < pdeData.unkGlobalIndex; p++)
					*(pdeData.unk1[p]) = pdeData.unk3[p];

				if (resMin > pdeData.threshold) {
					pdeData.maxDelta *= 1.5;
					selectSteps *= 2;
					tries++;//not acceptable, start over
				}
				else {
					break; // acceptable
				}
					
			} while (tries < 10);

			if (DBG_PDE_SOLVE && tries > 0) {
				pdeData.retriesCount += (double)100/steps;
				if (tries > pdeData.maxRetries)
					pdeData.maxRetries = tries;
			}

			if (tries == 10) {
				pdeData.score = HUGE_VAL;
				return false;
			}

			pdeData.maxDelta = initialDelta;
			pdeData.score += resMin / steps;
			
		}
	}
	return true;
}

Function** Function::SolvePDENL(int& maxSolutions, PDEEqn eqn, void* userdata, double delta, unsigned int selectSteps, unsigned int computeSteps, int precision, double rDelta, double valMin, double valMax) {

	double valSteps = Function::__Pow10(precision);
	Function** res = new Function*[maxSolutions];
	int k = 0;
	int tries = maxSolutions;

	double* scores = new double[tries];

	for (int i = 0; i < tries; i++){
		Function* testF = new Function(this, true);
		PDEData pdeData(*testF);
		pdeData.X = new double[_level];
		pdeData.s = new long long[_level];
		pdeData.UnitSteps = new double[_level];
		pdeData.valSteps = valSteps;
		pdeData.valMin = valMin;
		pdeData.valMax = valMax;
		pdeData.eqn = eqn;
		pdeData.n = _level;
		pdeData.userdata = userdata;
		pdeData.unkMaxSize = 2*_level * (_level+1);
		pdeData.unk1 = new double*[pdeData.unkMaxSize];
		pdeData.unk2 = new double[pdeData.unkMaxSize];
		pdeData.unk3 = new double[pdeData.unkMaxSize];
		pdeData.score = 0;
		pdeData.unkd = NULL;
		pdeData.maxDelta = delta;
		pdeData.rDelta = rDelta;
		pdeData.selectSteps = selectSteps;
		pdeData.computeSteps = computeSteps;

		SolvePDENL_MOLRecursive(0, pdeData);

		if (DBG_PDE_SOLVE)
			printf("solution found score: %lg\n", pdeData.score);
		if (k < maxSolutions){
			res[k] = testF;
			scores[k] = pdeData.score;
			k++;
		}
		else{ //no more place -- replace one with less score.
			for (int i = 0; i < k; i++)
				if (pdeData.score < scores[i]) {
				scores[i] = pdeData.score;
				Function* tf = res[i];
				res[i] = testF;
				testF = tf;
				}
			delete testF;
		}

	}

	for (int i = 0; i < k; i++){
		for (int j = 0; j < i; j++){
			if (scores[i] < scores[j]) {
				double tx = scores[i];
				Function* tf = res[i];
				scores[i] = scores[j];
				res[i] = res[j];
				scores[j] = tx;
				res[j] = tf;
			}
		}
	}
	delete[] scores;
	maxSolutions = k;
	return res;
}