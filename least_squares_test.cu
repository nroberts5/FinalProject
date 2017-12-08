// The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt
/*

    This is an example illustrating the use the general purpose non-linear 
    least squares optimization routines from the dlib C++ Library.

    This example program will demonstrate how these routines can be used for data fitting.
    In particular, we will generate a set of data and then use the least squares  
    routines to infer the parameters of the model which generated the data.
*/


#include "./dlib-19.7/dlib/optimization.h"
#include <iostream>
#include <vector>


using namespace std;
using namespace dlib;

// ----------------------------------------------------------------------------------------

typedef matrix<double,2,1> input_vector;
typedef matrix<double,3,1> parameter_vector;
typedef matrix<double,0,1> column_vector;

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

column_vector signal (const column_vector& unknowns, const column_vector& trs, const column_vector& tips, const column_vector& tes)
{


    const double beta = unknowns(0);
    const double T1__F = unknowns(1);
    const double T1__W = unknowns(2);
    const double rho__F = unknowns(3);
    const double rho__W = unknowns(4);
    const double R2s = unknowns(5);
    const double phi = unknowns(6);
    const double psi = unknowns(7);

    double TR;
    double alpha;
    double TE;
    double sR;
    double sI;

    std::vector<double> v;
    // compute signal model function and return the result
    for (int ACQNUM = 0; ACQNUM < trs.nr(); ACQNUM++)
    {
        for (int TENUM = 0; TENUM < tes.nr(); TENUM++)
        {
            TR = trs(ACQNUM);
            alpha = tips(ACQNUM);
            TE = tes(TENUM);
            sR = ((rho__W * sin(beta * alpha) * (1.0 - exp(-TR / T1__W)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * cos(2.0 * PI * fp[0] * TE) + a[1] * cos(2.0 * PI * fp[1] * TE) + a[2] * cos(2.0 * PI * fp[2] * TE) + a[3] * cos(2.0 * PI * fp[3] * TE) + a[4] * cos(2.0 * PI * fp[4] * TE) + a[5] * cos(2.0 * PI * fp[5] * TE)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * cos(TE * psi + phi) - rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * sin(2.0 * PI * fp[0] * TE) + a[1] * sin(2.0 * PI * fp[1] * TE) + a[2] * sin(2.0 * PI * fp[2] * TE) + a[3] * sin(2.0 * PI * fp[3] * TE) + a[4] * sin(2.0 * PI * fp[4] * TE) + a[5] * sin(2.0 * PI * fp[5] * TE)) * sin(TE * psi + phi) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
            sI = ((rho__W * sin(beta * alpha) * (1.0 - exp(-TR / T1__W)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * cos(2.0 * PI * fp[0] * TE) + a[1] * cos(2.0 * PI * fp[1] * TE) + a[2] * cos(2.0 * PI * fp[2] * TE) + a[3] * cos(2.0 * PI * fp[3] * TE) + a[4] * cos(2.0 * PI * fp[4] * TE) + a[5] * cos(2.0 * PI * fp[5] * TE)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * sin(TE * psi + phi) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * sin(2.0 * PI * fp[0] * TE) + a[1] * sin(2.0 * PI * fp[1] * TE) + a[2] * sin(2.0 * PI * fp[2] * TE) + a[3] * sin(2.0 * PI * fp[3] * TE) + a[4] * sin(2.0 * PI * fp[4] * TE) + a[5] * sin(2.0 * PI * fp[5] * TE)) * cos(TE * psi + phi) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
            v.push_back(sR);
            v.push_back(sI);
        }
    }
    
    matrix<double,0,1> sig = mat(v);
    return sig;
}

// We will use this function to generate data.  It represents a function of 2 variables
// and 3 parameters.   The least squares procedure will be used to infer the values of 
// the 3 parameters based on a set of input/output pairs.
double model (
    const input_vector& input,
    const parameter_vector& params
)
{
    const double p0 = params(0);
    const double p1 = params(1);
    const double p2 = params(2);

    const double i0 = input(0);
    const double i1 = input(1);

    const double temp = p0*i0 + p1*i1 + p2;

    return temp*temp;
}

// ----------------------------------------------------------------------------------------

// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double residual (
    const std::pair<input_vector, double>& data,
    const parameter_vector& params
)
{
    return model(data.first, params) - data.second;
}

double my_residual (
    const std::pair<input_vector, double>& data,
    const parameter_vector& params
)
{
    return model(data.first, params) - data.second;
}


// ----------------------------------------------------------------------------------------

// This function is the derivative of the residual() function with respect to the parameters.
parameter_vector residual_derivative (
    const std::pair<input_vector, double>& data,
    const parameter_vector& params
)
{
    parameter_vector der;

    const double p0 = params(0);
    const double p1 = params(1);
    const double p2 = params(2);

    const double i0 = data.first(0);
    const double i1 = data.first(1);

    const double temp = p0*i0 + p1*i1 + p2;

    der(0) = i0*2*temp;
    der(1) = i1*2*temp;
    der(2) = 2*temp;

    return der;
}

// ----------------------------------------------------------------------------------------

int main()
{
    // column_vector acq_params(3);
    // acq_params = 20e-3, 20*PI/180, 1.2e-3;

    int NPs = 3; // Number of TR/FlipAngle Pairs
    int NTEs = 6; // Number of Echoes
    int NACQS = NPs*NTEs;

    column_vector tes(NTEs); tes = 1.2e-3,3.2e-3,5.2e-3,7.2e-3,9.2e-3,11.2e-3;
    column_vector trs(NPs); trs = 5e-3,10e-3,15e-3;
    column_vector tips(NPs); tips = 6*PI/180,12*PI/180,80*PI/180;

    column_vector unknowns(8);
    unknowns = 1,312e-3,822e-3,50,950,30,0,0;

    cout<<signal(unknowns, trs, tips, tes);

    try
    {
        // randomly pick a set of parameters to use in this example
        const parameter_vector params = 10*randm(3,1);
        cout << "params: " << trans(params) << endl;


        // Now let's generate a bunch of input/output pairs according to our model.
        std::vector<std::pair<input_vector, double> > data_samples;
        input_vector input;
        for (int i = 0; i < 1000; ++i)
        {
            input = 10*randm(2,1);
            const double output = model(input, params);

            // save the pair
            data_samples.push_back(make_pair(input, output));
        }

        // Before we do anything, let's make sure that our derivative function defined above matches
        // the approximate derivative computed using central differences (via derivative()).  
        // If this value is big then it means we probably typed the derivative function incorrectly.
        cout << "derivative error: " << length(residual_derivative(data_samples[0], params) - 
                                               derivative(residual)(data_samples[0], params) ) << endl;





        // Now let's use the solve_least_squares_lm() routine to figure out what the
        // parameters are based on just the data_samples.
        parameter_vector x;
        x = 1;

        cout << "Use Levenberg-Marquardt" << endl;
        // Use the Levenberg-Marquardt method to determine the parameters which
        // minimize the sum of all squared residuals.
        solve_least_squares_lm(objective_delta_stop_strategy(1e-7).be_verbose(), 
                               residual,
                               residual_derivative,
                               data_samples,
                               x);

        // Now x contains the solution.  If everything worked it will be equal to params.
        cout << "inferred parameters: "<< trans(x) << endl;
        cout << "solution error:      "<< length(x - params) << endl;
        cout << endl;




        x = 1;
        cout << "Use Levenberg-Marquardt, approximate derivatives" << endl;
        // If we didn't create the residual_derivative function then we could
        // have used this method which numerically approximates the derivatives for you.
        solve_least_squares_lm(objective_delta_stop_strategy(1e-7).be_verbose(), 
                               residual,
                               derivative(residual),
                               data_samples,
                               x);

        // Now x contains the solution.  If everything worked it will be equal to params.
        cout << "inferred parameters: "<< trans(x) << endl;
        cout << "solution error:      "<< length(x - params) << endl;
        cout << endl;




        x = 1;
        cout << "Use Levenberg-Marquardt/quasi-newton hybrid" << endl;
        // This version of the solver uses a method which is appropriate for problems
        // where the residuals don't go to zero at the solution.  So in these cases
        // it may provide a better answer.
        solve_least_squares(objective_delta_stop_strategy(1e-7).be_verbose(), 
                            residual,
                            residual_derivative,
                            data_samples,
                            x);

        // Now x contains the solution.  If everything worked it will be equal to params.
        cout << "inferred parameters: "<< trans(x) << endl;
        cout << "solution error:      "<< length(x - params) << endl;

    }
    catch (std::exception& e)
    {
        cout << e.what() << endl;
    }
}

// Example output:
/*
params: 8.40188 3.94383 7.83099 

derivative error: 9.78267e-06
Use Levenberg-Marquardt
iteration: 0   objective: 2.14455e+10
iteration: 1   objective: 1.96248e+10
iteration: 2   objective: 1.39172e+10
iteration: 3   objective: 1.57036e+09
iteration: 4   objective: 2.66917e+07
iteration: 5   objective: 4741.9
iteration: 6   objective: 0.000238674
iteration: 7   objective: 7.8815e-19
iteration: 8   objective: 0
inferred parameters: 8.40188 3.94383 7.83099 

solution error:      0

Use Levenberg-Marquardt, approximate derivatives
iteration: 0   objective: 2.14455e+10
iteration: 1   objective: 1.96248e+10
iteration: 2   objective: 1.39172e+10
iteration: 3   objective: 1.57036e+09
iteration: 4   objective: 2.66917e+07
iteration: 5   objective: 4741.87
iteration: 6   objective: 0.000238701
iteration: 7   objective: 1.0571e-18
iteration: 8   objective: 4.12469e-22
inferred parameters: 8.40188 3.94383 7.83099 

solution error:      5.34754e-15

Use Levenberg-Marquardt/quasi-newton hybrid
iteration: 0   objective: 2.14455e+10
iteration: 1   objective: 1.96248e+10
iteration: 2   objective: 1.3917e+10
iteration: 3   objective: 1.5572e+09
iteration: 4   objective: 2.74139e+07
iteration: 5   objective: 5135.98
iteration: 6   objective: 0.000285539
iteration: 7   objective: 1.15441e-18
iteration: 8   objective: 3.38834e-23
inferred parameters: 8.40188 3.94383 7.83099 

solution error:      1.77636e-15
*/
