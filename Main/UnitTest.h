#pragma once

#ifndef π
#define π 3.14159265358979323846
#endif 
#define gauss(σ,μ,t) (exp(-(t-(μ))*(t-(μ))/(2*(σ)*(σ))) / ((σ)*sqrt(2 * π)))
#define gauss_cos(σ,μ,t,ω) (gauss(σ,μ,t) * cos(ω/(σ)*(t-(μ))) )
#define gauss_sin(σ,μ,t,ω) (gauss(σ,μ,t) * sin(ω/(σ)*(t-(μ))) )
#define gauss_cos2(σ,μ,t,ω) (gauss(σ,μ,t) * pow(cos(ω/(σ)*(t-(μ))),2) )
#define gauss_sin2(σ,μ,t,ω) (gauss(σ,μ,t) * pow(sin(ω/(σ)*(t-(μ))),2) )

class UnitTest
{
public:
	void FourierTest();
	void MultiArgTest();
	void SolvePDE();
	void OperatorsTest();
};

