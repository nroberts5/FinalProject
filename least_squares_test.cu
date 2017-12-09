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
#include <omp.h>
#include <chrono>
#include <random>
#include "tools.hpp"


using namespace std;
using namespace dlib;

// ----------------------------------------------------------------------------------------

typedef matrix<double,4,1> input_vector;
typedef matrix<double,8,1> parameter_vector;
typedef matrix<double,0,1> column_vector;

// ----------------------------------------------------------------------------------------

#define FIELD_STRENGTH 3.0
#define PI 3.141592654
#define STD_DEV 1.0

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);


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

double signal (const column_vector& acq_params, const column_vector& unknowns)
{
    const double beta = unknowns(0);
    const double T1__F = unknowns(1);
    const double T1__W = unknowns(2);
    const double rho__F = unknowns(3);
    const double rho__W = unknowns(4);
    const double R2s = unknowns(5);
    const double phi = unknowns(6);
    const double psi = unknowns(7);

    const double TR = acq_params(0);
    const double alpha = acq_params(1);
    const double TE = acq_params(2);
    const double real_part = acq_params(3);

    // compute signal model function and return the result
    double sR = ((rho__W * sin(beta * alpha) * (1.0 - exp(-TR / T1__W)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * cos(2.0 * PI * fp[0] * TE) + a[1] * cos(2.0 * PI * fp[1] * TE) + a[2] * cos(2.0 * PI * fp[2] * TE) + a[3] * cos(2.0 * PI * fp[3] * TE) + a[4] * cos(2.0 * PI * fp[4] * TE) + a[5] * cos(2.0 * PI * fp[5] * TE)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * cos(TE * psi + phi) - rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * sin(2.0 * PI * fp[0] * TE) + a[1] * sin(2.0 * PI * fp[1] * TE) + a[2] * sin(2.0 * PI * fp[2] * TE) + a[3] * sin(2.0 * PI * fp[3] * TE) + a[4] * sin(2.0 * PI * fp[4] * TE) + a[5] * sin(2.0 * PI * fp[5] * TE)) * sin(TE * psi + phi) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
    double sI = ((rho__W * sin(beta * alpha) * (1.0 - exp(-TR / T1__W)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * cos(2.0 * PI * fp[0] * TE) + a[1] * cos(2.0 * PI * fp[1] * TE) + a[2] * cos(2.0 * PI * fp[2] * TE) + a[3] * cos(2.0 * PI * fp[3] * TE) + a[4] * cos(2.0 * PI * fp[4] * TE) + a[5] * cos(2.0 * PI * fp[5] * TE)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * sin(TE * psi + phi) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * sin(2.0 * PI * fp[0] * TE) + a[1] * sin(2.0 * PI * fp[1] * TE) + a[2] * sin(2.0 * PI * fp[2] * TE) + a[3] * sin(2.0 * PI * fp[3] * TE) + a[4] * sin(2.0 * PI * fp[4] * TE) + a[5] * sin(2.0 * PI * fp[5] * TE)) * cos(TE * psi + phi) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
    
    return (real_part)*sR+(1-real_part)*sI;
}

double residual (const std::pair<column_vector, double>& data, const column_vector& params)
{
    return signal(data.first, params) - data.second;
}

std::vector<pair<input_vector,double>> signal_pairs(const int NACQS, const int NTES, const column_vector& trs, const column_vector& tips, const column_vector& tes, const column_vector& unknowns)
{
    std::vector<pair<input_vector,double>> data_samples;
    input_vector input;
    for (int nacq = 0; nacq < NACQS; nacq++)
    {
        for (int nte = 0; nte < NTES; nte++)
        {
            for (int real_part = 1; real_part >=0; real_part--)
            {
                input = trs(nacq), tips(nacq), tes(nte), real_part;
                const double sim_sig = signal(input,unknowns);//+distribution(generator); 
                data_samples.push_back(make_pair(input,sim_sig));
            }
        }
    }
    return data_samples;
}

std::vector<pair<input_vector,double>> add_noise2signal_pairs(const std::vector<pair<input_vector,double>> data_samples, const double stddev)
{
    normal_distribution<double> distribution(0.0,stddev);
    std::vector<pair<input_vector,double>> noisy_data_samples;
    for (auto ind_pair : data_samples)
    {
        const input_vector input = ind_pair.first;
        double sim_sig = ind_pair.second;
        noisy_data_samples.push_back(make_pair(input,sim_sig+distribution(generator)));
    }
    return noisy_data_samples;
}


column_vector NLSQ(std::vector<pair<input_vector,double>> noiseless_data_samples, const double std_dev_noise, const parameter_vector& unknowns)
{
    std::vector<pair<input_vector,double>> noisy_data_samples = add_noise2signal_pairs(noiseless_data_samples,std_dev_noise);
    parameter_vector x;
    x=1, 300e-3, 800e-3, 5000, 5000, 0, 0, 0;
    // cout << "Use Levenberg-Marquardt, approximate derivatives" << endl;
    solve_least_squares_lm(objective_delta_stop_strategy(1e-10)/*.be_verbose()*/, 
                           residual,
                           derivative(residual),
                           noisy_data_samples,
                           x);
    return x;
    // cout << "inferred parameters: "<< trans(x) << endl;
    // cout << "solution error:      "<< length(x - unknowns) << endl;
    // cout << endl;
}

int main(int argc, char const *argv[])
{
    int NUMTHREADS;
    int NSIMS;
    int NUMTRIALS;

    if (argc<3)
    {
        cout<<"Usage ./main NUMTHREADS NSIMS NUMTRIALS\n";
        return 1;
    }
    else
    {
        NUMTHREADS = atoi(argv[1]);
        NSIMS = atoi(argv[2]);
        NUMTRIALS = atoi(argv[3]);

    }

    int NACQS = 3; // Number of TR/FlipAngle Pairs
    int NTES = 6; // Number of Echoes

    column_vector trs(NACQS); trs = 5e-3,10e-3,15e-3;
    column_vector tips(NACQS); tips = 6*PI/180,12*PI/180,80*PI/180;
    column_vector tes(NTES); tes = 1.2e-3,3.2e-3,5.2e-3,7.2e-3,9.2e-3,11.2e-3;

    double beta=1.0, T1__F=312e-3, T1__W=822e-3, rho__F=50, rho__W=9050, R2s=30, phi=0, psi=0;
    parameter_vector unknowns; unknowns = beta, T1__F, T1__W, rho__F, rho__W, R2s, phi, psi;

    std::vector<pair<input_vector,double>> noiseless_data_samples = signal_pairs(NACQS, NTES, trs, tips, tes, unknowns);
    
    struct stopwatch sw;

    double betas[NSIMS]{};

    std::vector<double> timings;
    for (int trial = 0; trial < NUMTRIALS; trial++)
    {
        sw.click();
        #pragma omp parallel num_threads(NUMTHREADS)
        {
            #pragma omp for nowait
            for (int i = 0; i < NSIMS; i++)
            {
                column_vector x = NLSQ(noiseless_data_samples, 1.0, unknowns);
                betas[i] = x(0);
            }
        }
        sw.click();
        timings.push_back(sw.check());
    }
    column_vector t = mat(timings);
    double minTIME = dlib::min(t);
    cout<<NSIMS<<", "<<NUMTHREADS<<", "<<minTIME/1000.0<<"s"<<", <"<<NUMTRIALS<<" trial(s)>"<<endl;

    return 0;
}
