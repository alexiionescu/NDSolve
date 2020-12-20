
#include <stdio.h>
#include <stdlib.h>
#include "Plot.h"
#include "Dislin\dislin_d.h"


#define MAX_PLOT_POINTS	2000
#define MAX_COLORS	6
#define MAX_LEGEND_LINE	100
void __Plot(int N, Function ** y, const char** titleLabel, const char* xAxisLabel, const char* yAxisLabel, const char** legendLabels, const char* cdev, bool bParametric)
{
	if (N < 1 || !y || (bParametric && N % 2 != 0))
		return;

	//find min and max time
	double MinT = HUGE_VAL, MaxT = -HUGE_VAL;
	for (int k = 0; k<N; k++){
		if (MinT >(*y[k]).GetT0()) MinT = (*y[k]).GetT0();
		if (MaxT < (*y[k]).GetT())  MaxT = (*y[k]).GetT();
	}

	long long n = (*y[0]).GetSize();
	double* times = bParametric ? NULL : (double*)malloc(MAX_PLOT_POINTS*sizeof(double));
	double* yval = (double*)malloc(N*MAX_PLOT_POINTS*sizeof(double));

	double MinY = HUGE_VAL, MaxY = -HUGE_VAL;
	double MinX, MaxX;
	if (bParametric) {
		MinX = HUGE_VAL;
		MaxX = -HUGE_VAL;
	}
	else {
		MinX = MinT;
		MaxX = MaxT;
	}

	for (int gn = 0; gn < MAX_PLOT_POINTS; gn++)
	{
		double t = MinT + gn * (MaxT - MinT) / MAX_PLOT_POINTS;
		if (!bParametric) times[gn] = t;
		for (int k = 0; k<N; k++){
			double x = (*y[k])[t];
			yval[k*MAX_PLOT_POINTS + gn] = x;
			if (!bParametric || k % 2 == 1)
			{
				if (MinY > x) MinY = x;
				if (MaxY < x) MaxY = x;
			}
			else //for parametric recalculate MinX and MaxX based on functions 0,2,4,...
			{
				if (MinX > x) MinX = x;
				if (MaxX < x) MaxX = x;
			}
		}
	}

	
	if (titleLabel && titleLabel[0])
		setfil(titleLabel[0]);
	metafl(cdev); //can use BMP, GIF,etc
	setpag("DA4L");
	if (0 == strcmp(cdev, "XWIN")) {
		scrmod("AUTO");
		int nw, nh;
		getscr(&nw, &nh);
		window(nw/6, nh/6, 2 * nw / 3, 3* nh / 4);
	} else {
		scrmod("REVERSE");
		winsiz(2969, 2099);
	}
	disini ();
	  pagera ();
	  hwfont ();
	  
	  axspos (350, 1900);
	  axslen (2440, 1650);
	  labdig (2, "x");
	  labdig (2, "y");
      ticks  (10, "xy");
	  
	  texmod("ON"); // when set you can put text in TeX form for example $\\frac{1}{\\sqrt{x^2+y^2}}$
							
	  height(22);
	  hname(22);
	  if(titleLabel){
		  for(int i=0; i<4 && titleLabel[i];i++)
		  	  titlin(titleLabel[i],i+1);
	  }
		  
	  char* legbuf = NULL;
	  if(legendLabels)
	  {
		  legbuf = new char[(bParametric ? N/2 : N)*MAX_LEGEND_LINE];
		  legini(legbuf,bParametric ? N/2 : N,MAX_LEGEND_LINE);
	  }

	  if(xAxisLabel) name(xAxisLabel, "x");
	  if(yAxisLabel) name(yAxisLabel, "y");	  
	   
	  double colors[MAX_COLORS][3]={{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1}};

	   graf   (MinX, MaxX, MinX, (MaxX-MinX)/10, 
		       MinY, MaxY, MinY, (MaxY-MinY)/10);
	   

	   if(titleLabel) title();
	   for(int k=0;k<N;k++)
	   {
		   setrgb(colors[k%MAX_COLORS][0],colors[k%MAX_COLORS][1],colors[k%MAX_COLORS][2]);
		   if(bParametric)
		   {
			   curve  (yval+k*MAX_PLOT_POINTS, yval+MAX_PLOT_POINTS*(k+1), MAX_PLOT_POINTS);
			   if(legendLabels){ 
				   if(legendLabels[k/2]) 
					   leglin(legbuf,legendLabels[k/2],k/2+1);
				   else legendLabels = NULL;
			   }

			   k++;
		   } else {
			   curve  (times, yval+k*MAX_PLOT_POINTS, MAX_PLOT_POINTS);
			   if(legendLabels) {
				   if(legendLabels[k]) 
					   leglin(legbuf,legendLabels[k],k+1);
				   else legendLabels = NULL;
			   }
		   }
	   }

	  color  ("fore");
	  if(legendLabels) {
		  height(24);
		  legend(legbuf,4);
	  }
	  if(legbuf) free(legbuf);
	  dash   ();
	  xaxgit ();
	  yaxgit ();
   disfin ();

    free(times);
	free(yval);
}

void Plot(int N,Function ** y,const char** title, const char* xAxisLabel,const char* yAxisLabel,const char** legendLabels,const char* cdev)
{
	__Plot(N,y,title,xAxisLabel,yAxisLabel,legendLabels,cdev,false);
}
void ParametricPlot(int N,Function ** y,const char** title, const char* xAxisLabel,const char* yAxisLabel,const char** legendLabels,const char* cdev)
{
	__Plot(N,y,title,xAxisLabel,yAxisLabel,legendLabels,cdev,true);
}