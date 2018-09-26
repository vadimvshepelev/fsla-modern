#include<stdio.h>

int hello(void) {
    printf("Hello, world!\n");
    return 6;
}

double fL(double gamma, double p, double roL, double vL, double pL) {	
	double f = 0.;
	if(p>pL) {
		double AL = 2./(gamma+1)/roL;
		double BL = (gamma-1.)/(gamma+1.)*pL;
		f = (p-pL) * sqrt(AL/(p+BL));
		return f;
	} else {
		double cL = sqrt(gamma*pL/roL);
		f = 2.*cL/(gamma-1.) * ( (pow(p/pL, (gamma-1.)/2./gamma)) - 1. );
		return f;	
	}
}

double fR(double gamma, double p, double roR, double vR, double pR) {
	double f = 0.;
	if(p>pR) {
		double AR = 2./(gamma+1)/roR;
		double BR = (gamma-1.)/(gamma+1.)*pR;
		f = (p-pR) * sqrt(AR/(p+BR));
		return f;
	} else {
		double cR = sqrt(gamma*pR/roR);
		f = 2.*cR/(gamma-1.) * ( (pow(p/pR, (gamma-1.)/2./gamma)) - 1. );
		return f;
	}
}

