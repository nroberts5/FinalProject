// The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt
/*

    This is an example illustrating the use the general purpose non-linear 
    optimization routines from the dlib C++ Library.

    The library provides implementations of the conjugate gradient,  BFGS,
    L-BFGS, and BOBYQA optimization algorithms.  These algorithms allow you to
    find the minimum of a function of many input variables.  This example walks
    though a few of the ways you might put these routines to use.

*/


#include "./dlib-19.7/dlib/optimization.h"
#include <iostream>


using namespace std;
using namespace dlib;

// ----------------------------------------------------------------------------------------

// In dlib, the general purpose solvers optimize functions that take a column
// vector as input and return a double.  So here we make a typedef for a
// variable length column vector of doubles.  This is the type we will use to
// represent the input to our objective functions which we will be minimizing.
typedef matrix<double,0,1> column_vector;

// ----------------------------------------------------------------------------------------
// Below we create a few functions.  When you get down into main() you will see that
// we can use the optimization algorithms to find the minimums of these functions.
// ----------------------------------------------------------------------------------------

#define FIELD_STRENGTH 3.0
#define PI 3.141592654
float ppm[6] = {0.6, -0.5, -1.95, -2.60, -3.40, -3.80};

float* fp_maker(const float* ppm)
{
    float* fp = new float[6]; 
    memcpy(fp,ppm,6*sizeof(float)); 
    for(int i = 0; i < 6; i++)
    {
        fp[i] *= (FIELD_STRENGTH*42.577);
    }
    return fp;
}

float* fp = fp_maker(ppm);
float a[6] = {4.7/100, 3.9/ 100, 0.6/ 100, 12.0/ 100, 70.0/ 100, 8.8/ 100};



float* signal_equation(const float TR, const float alpha, const float TE, const float beta, const float T1__F, const float T1__W, const float rho__F, const float rho__W, const float R2s, const float phi, const float psi)
{   
    float sR = ((rho__W * sin(beta * alpha) * (1.0 - exp(-TR / T1__W)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * cos(2.0 * PI * fp[0] * TE) + a[1] * cos(2.0 * PI * fp[1] * TE) + a[2] * cos(2.0 * PI * fp[2] * TE) + a[3] * cos(2.0 * PI * fp[3] * TE) + a[4] * cos(2.0 * PI * fp[4] * TE) + a[5] * cos(2.0 * PI * fp[5] * TE)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * cos(TE * psi + phi) - rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * sin(2.0 * PI * fp[0] * TE) + a[1] * sin(2.0 * PI * fp[1] * TE) + a[2] * sin(2.0 * PI * fp[2] * TE) + a[3] * sin(2.0 * PI * fp[3] * TE) + a[4] * sin(2.0 * PI * fp[4] * TE) + a[5] * sin(2.0 * PI * fp[5] * TE)) * sin(TE * psi + phi) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
    float sI = ((rho__W * sin(beta * alpha) * (1.0 - exp(-TR / T1__W)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * cos(2.0 * PI * fp[0] * TE) + a[1] * cos(2.0 * PI * fp[1] * TE) + a[2] * cos(2.0 * PI * fp[2] * TE) + a[3] * cos(2.0 * PI * fp[3] * TE) + a[4] * cos(2.0 * PI * fp[4] * TE) + a[5] * cos(2.0 * PI * fp[5] * TE)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * sin(TE * psi + phi) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * sin(2.0 * PI * fp[0] * TE) + a[1] * sin(2.0 * PI * fp[1] * TE) + a[2] * sin(2.0 * PI * fp[2] * TE) + a[3] * sin(2.0 * PI * fp[3] * TE) + a[4] * sin(2.0 * PI * fp[4] * TE) + a[5] * sin(2.0 * PI * fp[5] * TE)) * cos(TE * psi + phi) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
    
    float* sig = new float[2];
    sig[0] = sR;
    sig[1] = sI;
    return sig;
}

float* signal (const column_vector& params, const column_vector& acq)
{
    const double TR = acq(0);
    const double alpha = acq(1);
    const double TE = acq(2);

    const double beta = params(0);
    const double T1__F = params(1);
    const double T1__W = params(2);
    const double rho__F = params(3);
    const double rho__W = params(4);
    const double R2s = params(5);
    const double phi = params(6);
    const double psi = params(7);


    // compute signal model function and return the result
    float sR = ((rho__W * sin(beta * alpha) * (1.0 - exp(-TR / T1__W)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * cos(2.0 * PI * fp[0] * TE) + a[1] * cos(2.0 * PI * fp[1] * TE) + a[2] * cos(2.0 * PI * fp[2] * TE) + a[3] * cos(2.0 * PI * fp[3] * TE) + a[4] * cos(2.0 * PI * fp[4] * TE) + a[5] * cos(2.0 * PI * fp[5] * TE)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * cos(TE * psi + phi) - rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * sin(2.0 * PI * fp[0] * TE) + a[1] * sin(2.0 * PI * fp[1] * TE) + a[2] * sin(2.0 * PI * fp[2] * TE) + a[3] * sin(2.0 * PI * fp[3] * TE) + a[4] * sin(2.0 * PI * fp[4] * TE) + a[5] * sin(2.0 * PI * fp[5] * TE)) * sin(TE * psi + phi) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
    float sI = ((rho__W * sin(beta * alpha) * (1.0 - exp(-TR / T1__W)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * cos(2.0 * PI * fp[0] * TE) + a[1] * cos(2.0 * PI * fp[1] * TE) + a[2] * cos(2.0 * PI * fp[2] * TE) + a[3] * cos(2.0 * PI * fp[3] * TE) + a[4] * cos(2.0 * PI * fp[4] * TE) + a[5] * cos(2.0 * PI * fp[5] * TE)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * sin(TE * psi + phi) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * sin(2.0 * PI * fp[0] * TE) + a[1] * sin(2.0 * PI * fp[1] * TE) + a[2] * sin(2.0 * PI * fp[2] * TE) + a[3] * sin(2.0 * PI * fp[3] * TE) + a[4] * sin(2.0 * PI * fp[4] * TE) + a[5] * sin(2.0 * PI * fp[5] * TE)) * cos(TE * psi + phi) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
    
    float* sig = new float[2];
    sig[0] = sR;
    sig[1] = sI;
    return sig;
}

// ----------------------------------------------------------------------------------------

class test_function
{
public:

    test_function (const column_vector& input)
    {
        target = input;
    }

    double operator() ( const column_vector& arg) const
    {
        // return the mean squared error between the target vector and the input vector
        return mean(squared(target-arg));
    }

private:
    column_vector target;
};

// ----------------------------------------------------------------------------------------


int main(int argc, char const *argv[])
{
    /* code */
    column_vector target(4);
    target = 1,2,3,4;
    struct test_function tf = test_function(target);
    target = 2,3,4,5;
    cout<<tf(target)<<endl;

    return 0;
}

