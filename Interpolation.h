#ifndef INTERPOLATION_H
#define INTERPOLATION_H


class Interpolation
{
    public:
        Interpolation();
        virtual ~Interpolation();
        /**
        @param x[]  all points
        @param y[]  corresponding f(x)
        @param a    point to be interpolated
        @param n size of the input array
        @return interpolated point
        */
        double NewtInt(double x[], double y[], double a, int n);
        /**
        @param x[]  all points
        @param y[]  corresponding f(x)
        @param a    point to be interpolated
        @param n size of the input array
        @return interpolated point
        */
        double Spline(double x[], double y[], double a, int n);
    protected:

    private:
    	///////////////Spline Interpolation/////////////////
    	
		/**
    	@param x[] array of interval limits
    	@param a   point to be interpolated
    	@param n   number of elements in array x[]
    	
    	Finds the interval of the interpolated point to calculate the interpolant
    	*/
    	double findingInterval(double x[], double a, int n);
    	
    	
    	
    	//////////////Newton Interpolation//////////////////////
    	
		/**
    	@param a point to be interpolated
    	@param fdd[] array of finite divided differences acting as coefficients for the equation
    	@param n number of points given for training data
    	@param x input array containing the training data
    	
		Interpolates the point a given the finite divided differences according to the following equation
		fn(x) = b0 + b1(x - x0) +···+ bn(x - x0)(x - x1)···(x - xn-1)
		whereas b0, b1..  are the finite divided differences
    	*/
    	double evaluatingNewtonPoint(double a, double fdd[], int n, double x[]);
    	/**
    	@param x[] array of input values 
    	@param y[] array of output values for initializing fdds
    	@param n size of the array
    	*/
    	double* calculateFdds(double x[], double y[], int n);
};

#endif // INTERPOLATION_H
