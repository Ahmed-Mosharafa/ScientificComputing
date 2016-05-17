#include "Interpolation.h"
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;

Interpolation::Interpolation()
{
    //ctor
}

Interpolation::~Interpolation()
{
    //dtor
}
////////////// Spline Interpolation ////////////////////////

double Interpolation::findingInterval(double x[], double a, int n) {
		for (int i = n-1; i > 0; i--){
			if (a < x[i] && a > x[i-1])	{
				return i-1;	
			}
		} 
}

double Interpolation::Spline(double x[], double y[], double a, int n){
	int interval, i;
	double h[n], alpha[n], l[n], u[n], z[n], c[n], b[n], d[n], interpolant,aj,bj,cj,dj;
	//setting length of the intervals 
	for (i=0; i<n-1; i++){	
		h[i] = x[i+1]-x[i];  	
	}
	// 
	for (i=0; i<n-1; i++){
		alpha[i] = ((3/h[i]) * (y[i+1] - y[i]) )- ( (3/h[i-1]) * (y[i] - y[i-1]) );
	}
	//Tridiagonal matrix formation
	//z[0], z[n] = 0
	l[0] = 1;
 		u[0] = 0;
	z[0] = 0;
	for (i = 1; i<n-1; i++){
		l[i] = (2 * (x[i+1] - x[i-1])) - (h[i-1] * u[i-1]);
		u[i] = h[i] / l[i]; 
		z[i] = (alpha[i] - (h[i-1] * z[i-1]) ) / l[i]; //second derivative of f(x[i])
	}
	l[n] = 1;
	z[n] = 0;
	c[n] = 0;
	for (i=n-1; i>=0; i--){
		c[i] = z[i] - (u[i] * c[i+1]) ;
		b[i] = ( (y[i+1] - y[i]) / (h[i]) ) -  ((h[i]/3) * (c[i+1] + (2*c[i])));
		d[i] = (c[i+1] - c[i]) / (3*h[i]);
	}
	sort(x,x+n);
	if (a < x[0] || a > x[n-1]){
		cout << "Not in range";
		return 0;
	}
	else{
		interval = findingInterval(x, a, n);
		//calculated coefficients of Sj(x) = aj + bj(x-xj) + cj(x-xj)^2 + dj(x-xj) ^3
		aj = y[interval];
		bj = b[interval];
		cj = c[interval];
		dj = d[interval];
		interpolant = aj + bj + pow(  cj * (a - x[interval]) , 2) + pow( dj * (a - x[interval]) , 3);
		cout << "\naj = " << aj << "\nbj = " << bj << "\ncj = " << cj << "\ndj = " << dj << "\n" << "x[interval] = " << x[interval] << "\n	" << "a: " << a << "\n";
		return interpolant;
	}
}

double Interpolation::evaluatingNewtonPoint(double a, double fdd[], int n, double x[]){
	double sum, mult;
	for(int i=n-1;i>=0;i--)
    {
        mult = 1;
        for(int j=0;j<i;j++){
            mult *= (a-x[j]);
        }
        mult *= fdd[i];
        sum  += mult;
    cout << "Iteration " << (n-1)-i << ": " << sum << "\n";
    }
    return sum;
}
double* Interpolation::calculateFdds(double x[], double y[], int n){
    double fdd[n];
    //initializing FDDs with Ys 
    for (int i=0; i<n; i++){
        fdd[i] = y[i];
    }
    //calculate Fdds 
	for (int j=0; j<n-1; j++){
        for (int i=n-1; i>j;i--){	
            	fdd[i] = ( fdd[i]-fdd[i-1] )/ ( x[i]-x[i-j-1]) ;
        }
    }
    return fdd;
}
double Interpolation::NewtInt(double x[], double y[], double a, int n)
{
    double *fdd;   //initializng an array of Fdds 
	double interpolant;
	int i,j;
	
    //calculating FDDs for each iteration 
	fdd = calculateFdds(x, y , n);
	
    //evaluating fn(x) while calculating coefficients 
    interpolant = evaluatingNewtonPoint(a, fdd, n, x);
    return interpolant;
}

/**
	@param a[] array to be printed
	@param n   size of the array
	helper function for printing the array
*/
void printArray(double a[], int n){
	for (int i = 0; i<n; i++){
		cout << a[i] << " ";
	}
	cout << "\n";
}
int main(){
    double newton_test[]   = {1,4,6,5,3,1.5,2.5}; //3.5
    double newton_values[] = {0.0,1.3862944,1.7917595,1.6094379,1.0986123,0.4054641,0.9162907}; //1.2527630
    double spline_test[]   = {3,4.5,7};
	double spline_values[] = {2.5,1,2.5};
	double spline_point    = 5;	
    double newton_point    = 3.5;
	Interpolation interpolate ;  
	cout << "Test Data of Newton Interpolation\n"; 
	printArray(newton_test, 7);
	printArray(newton_values, 7);
	cout << "Newton Interpolation of point " << newton_point << " is: " << interpolate.NewtInt(newton_test, newton_values, newton_point ,7) << "\n\n";
	cout << "Test Data of Cubic Spline\n"; 
	printArray(spline_test, 3);
	printArray(spline_values, 3);
	cout << "\nSpline Interpolation of point " << spline_point << " is: " << interpolate.Spline(spline_test, spline_values, spline_point ,3);
	return 0;
}
