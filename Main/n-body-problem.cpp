#include <stdio.h>
#include <stdlib.h>

#include "Function.h"
using namespace FunctionalMath;

#include "Plot.h"
#include "n-body-problem.h"


/*defines for n-body system equations of motions k=0..n-1 
Vx0' = f[0](t,x0,Vx0,y0,Vy0,....) // this is acceleration on direction x for first body -> Obtained from Ax = sum(Fx) / m newton's law for direction x of the first body
x0' = f[1] -> Vx0 // this is speed on direction x law .. all F[2k+1] Function * just return the 
Vy0' = f[2](t,x0,Vx0,y0,Vy0,....)
y0' = f[3] -> Vy0
Vz0' = f[4](t,x0,Vx0,y0,Vy0,....)
z0' = f[5] -> Vz0
...
Vxk' = f[6k](t,x1,Vx1,y1,Vy1,....)
xk' = f[6k+1] -> Vxk 
Vyk' = f[6k+2](t,x1,Vx1,y1,Vy1,....)
yk' = f[6k+3] -> Vyk
Vzk' = f[6k+4](t,x1,Vx1,y1,Vy1,....)
zk' = f[6k+5] -> Vzk
...
*/
#define Vx(k)	F(6*k)
#define x(k)	F(6*k+1)
#define Vy(k)	F(6*k+2)
#define y(k)	F(6*k+3)
#define Vz(k)	F(6*k+4)
#define z(k)	F(6*k+5)

static const double days=10*365.25; //days to compute


#define SecPerDay 86400
#define G 4.98137E-10   //gravitational constant in km^3/(kg days^2)




double Distance(int i,int j,FunctionSystem& F)
{
	return sqrt(pow(x(i)-x(j),2)+
				pow(y(i)-y(j),2)+
				pow(z(i)-z(j),2));
}

#define FuncAcc(i,Q) Acc##Q##i
static double *Gm;
static int n;
//this is equation motion a=F/m for axis Q and body i 
//use if(#@Q=='x' && i==0) return a+2; to modify
#define Acc(i,Q) double Acc##Q##i(double t, FunctionSystem& F,void*){\
	double a=0;\
	for(int j=0; j<n;j++)\
	{\
		if(i!=j) {\
		a += Gm[j]*(##Q(j)-##Q(i))/pow(Distance(i,j,F),3);\
		}\
	}\
	return a;\
}


#define FuncSpeed(i,Q) Speed##Q##i
#define Speed(i,Q) double Speed##Q##i(double t,FunctionSystem& F,void*) { return V##Q(i);}

#define NameEqn(i) FuncAcc(i,x),FuncSpeed(i,x),FuncAcc(i,y),FuncSpeed(i,y),FuncAcc(i,z),FuncSpeed(i,z)
#define Eqn(i)     Acc(i,x) Speed(i,x) Acc(i,y) Speed(i,y) Acc(i,z) Speed(i,z)

Eqn(0)
Eqn(1)
Eqn(2)
Eqn(3)
Eqn(4)
Eqn(5)

void NBodyProblem::NewObjectSolarSystem()
{
	double localGm[] = {
		G* 1.989E30,  //Sun
		G* 5.972E24,  //Earth
		G* 3E14,  //Comet
	};
	Gm = localGm;
	const int localN = (sizeof(localGm) / sizeof(localGm[0]));
	n = localN;

	double initialValues[6 * localN] = {
		/*Vxk[0] km/day,	xk[0] km,		Vyk[0] km/day,	yk[0] km,		Vzk[0] km/day,	zk[0] km*/
		0,					0,				0,				0,				0,				0,	 //Sun is center
		0,					148E6,			30*SecPerDay,	0,				0,				SecPerDay,	 //Sun is center
		-1*SecPerDay,		1E9,			3*SecPerDay,	0,			    1*SecPerDay,	0,	 //Comet
	};
	FunctionSystem::FunctionDef odeFunc[6 * localN] = {
		NameEqn(0),
		NameEqn(1),
		NameEqn(2)
	};

	FunctionSystem fs(n*6,days,initialValues,0,2);
	int tick = os_ticks();
	
	fs.SolveODENL(odeFunc);
	
	int tickn = os_ticks();
	printf("Solve Done in %d ms\n",tickn-tick);
	tick = tickn;
	
	int ref=0,obs=1,nobs=2;
	const char* titles[] = {"Comet","",NULL};

	
	Function * Xro=fs.CalculateFunction([ref,obs](double t, FunctionSystem& F) {return (x(obs)-x(ref))/1E6;});
	Function * Yro=fs.CalculateFunction([ref,obs](double t, FunctionSystem& F) {return (y(obs)-y(ref))/1E6;});
	Function * Zro=fs.CalculateFunction([ref,obs](double t, FunctionSystem& F) {return (z(obs)-z(ref))/1E6;});
	Function * Xrn=fs.CalculateFunction([ref,nobs](double t, FunctionSystem& F) {return (x(nobs)-x(ref))/1E6;});
	Function * Yrn=fs.CalculateFunction([ref,nobs](double t, FunctionSystem& F) {return (y(nobs)-y(ref))/1E6;});
	Function * Zrn=fs.CalculateFunction([ref,nobs](double t, FunctionSystem& F) {return (z(nobs)-z(ref))/1E6;});
	Function *  funcToPlotp[]={Xro,Yro,Xrn,Yrn};


	const char* plegends[] = {"XY earth","XY New Object",NULL};
	ParametricPlot(sizeof(funcToPlotp)/sizeof(funcToPlotp[0]),funcToPlotp,titles,"mil. km","mil km",plegends);
	delete Xro;delete Yro;delete Zro;delete Xrn;delete Yrn;delete Zrn;
	
	
	Function * Rno=fs.CalculateFunction([nobs,obs](double t, FunctionSystem& F) {return Distance(nobs,obs,F)/1E6;});
	Function * Rrn=fs.CalculateFunction([ref,nobs](double t, FunctionSystem& F) {return Distance(ref,nobs,F)/1E6;});
	Function *  funcToPlot[]={Rno,Rrn};
	const char* legends[] = {"Dist. Earth-Comet","Dist. Sun-Comet",NULL};

	AnalyzesResult a{};
	for(int i=0; legends[i];i++)
	{
		printf("\n------ %s -------\n",legends[i]);
		funcToPlot[i]->Analyze(a,stdout);
	}
	Plot(sizeof(funcToPlot)/sizeof(funcToPlot[0]),funcToPlot,titles,"days","",legends);
	delete Rno;
	delete Rrn;
}