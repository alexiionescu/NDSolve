#pragma once
#include "Function.h"
using namespace FunctionalMath;
//Plot 2D
//can provide text in TeX form for example: $\\frac{1}{\\sqrt{x^2+y^2}}$
//cdev can be BMP,GIF, etc...
//title and legendLabels are null terminated arrays
void Plot(int N,Function **y,const char** title=NULL, const char* xAxisLabel=NULL,const char* yAxisLabel=NULL,const char** legendLabels=NULL,const char* cdev="XWIN");


//Parametrix Plot 2D y[2k] on X-axis and Y[2k+1] on Y-axis
//parameters are same as Plot 2D
void ParametricPlot(int N,Function **y,const char** title=NULL, const char* xAxisLabel=NULL,const char* yAxisLabel=NULL,const char** legendLabels=NULL,const char* cdev="XWIN");
