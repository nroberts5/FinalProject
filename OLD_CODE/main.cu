#include <iostream>
#include <cstring>
#include <cmath>
// #include "H5Cpp.h"
#include <vector>
#include <random>
#include <chrono>
#include <thrust/transform.h>
#include <omp.h>
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/normal_distribution.h>
#include <thrust/generate.h>
#include <thrust/host_vector.h>
#include <thrust/execution_policy.h>
#include <cstdlib>
#include "./dlib-19.7/dlib/optimization.h"
#include "tools.hpp"


using namespace std;
// using namespace dlib;


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



float* signal_equation(const float TR, const float alpha, const float TE, const float beta, const float T1__F, const float T1__W, const float rho__F, const float rho__W, const float R2s, const float phi, const float psi){	
	float sR = ((rho__W * sin(beta * alpha) * (1.0 - exp(-TR / T1__W)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * cos(2.0 * PI * fp[0] * TE) + a[1] * cos(2.0 * PI * fp[1] * TE) + a[2] * cos(2.0 * PI * fp[2] * TE) + a[3] * cos(2.0 * PI * fp[3] * TE) + a[4] * cos(2.0 * PI * fp[4] * TE) + a[5] * cos(2.0 * PI * fp[5] * TE)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * cos(TE * psi + phi) - rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * sin(2.0 * PI * fp[0] * TE) + a[1] * sin(2.0 * PI * fp[1] * TE) + a[2] * sin(2.0 * PI * fp[2] * TE) + a[3] * sin(2.0 * PI * fp[3] * TE) + a[4] * sin(2.0 * PI * fp[4] * TE) + a[5] * sin(2.0 * PI * fp[5] * TE)) * sin(TE * psi + phi) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
	float sI = ((rho__W * sin(beta * alpha) * (1.0 - exp(-TR / T1__W)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * cos(2.0 * PI * fp[0] * TE) + a[1] * cos(2.0 * PI * fp[1] * TE) + a[2] * cos(2.0 * PI * fp[2] * TE) + a[3] * cos(2.0 * PI * fp[3] * TE) + a[4] * cos(2.0 * PI * fp[4] * TE) + a[5] * cos(2.0 * PI * fp[5] * TE)) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * sin(TE * psi + phi) + rho__F * sin(beta * alpha) * (1.0 - exp(-TR / T1__F)) * (a[0] * sin(2.0 * PI * fp[0] * TE) + a[1] * sin(2.0 * PI * fp[1] * TE) + a[2] * sin(2.0 * PI * fp[2] * TE) + a[3] * sin(2.0 * PI * fp[3] * TE) + a[4] * sin(2.0 * PI * fp[4] * TE) + a[5] * sin(2.0 * PI * fp[5] * TE)) * cos(TE * psi + phi) / (1.0 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
	
	float* sig = new float[2];
	sig[0] = sR;
	sig[1] = sI;
	return sig;
}

float* noiseless_signal_vector(const int N, const float* TR, const float* alpha, const float* TE, const float beta, const float T1__F, const float T1__W, const float rho__F, const float rho__W, const float R2s, const float phi, const float psi, const int NUMTHREADS)
{
	float* vect_out = new float[2*N]{};
	float* sig;

	#pragma omp parallel shared(vect_out) private(sig) num_threads(NUMTHREADS)
	{
		#pragma omp for nowait
		for(int n=0; n<N; n+=2)
		{
			sig = signal_equation(TR[n], alpha[n], TE[n], beta, T1__F, T1__W, rho__F, rho__W, R2s, phi, psi);
			vect_out[n] = sig[0];
			vect_out[n+1] = sig[1];
		}

	}
	return vect_out;
}

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

float* add_noise2vector(const float* vect_in, const int N, const float stddev, const int NUMTHREADS)
{  
	float* vect_out = new float[N]{}; 
	normal_distribution<float> distribution(0.0,1.0);
	#pragma omp parallel shared(vect, stddev) private(seed) num_threads(NUMTHREADS)
	{
		#pragma omp for nowait
		for (int n = 0; n < N; n++)
		{
			vect_out[n]=vect_in[n]+distribution(generator);
		}

	}
	return vect_out;
}

float* tile_vector(const float* vect, const int n, const int numREPS, const int NUMTHREADS)
{
	float* vect_out = new float[numREPS*n];
	
	#pragma omp parallel shared(vect_out, numREPS, n) num_threads(NUMTHREADS)
	{
		#pragma omp for nowait
		for (int i = 0; i < numREPS; i++)
		{
			memcpy(&vect_out[i*n], vect, n*sizeof(*vect));
		}
	}

	return vect_out;
}


int main(int argc, char const *argv[])
{
	int NUMTHREADS;
	int NSIMS;

	if (argc<3)
	{
		cout<<"Usage ./main NUMTHREADS NSIMS\n";
		return 1;
	}
	else
	{
		NUMTHREADS = atoi(argv[1]);
		NSIMS = atoi(argv[2]);
	}

	int NPs = 3; // Number of TR/FlipAngle Pairs
	int NTEs = 6; // Number of Echoes
	int NACQS = NPs*NTEs;

	float tes[NTEs] = {1.2e-3,3.2e-3,5.2e-3,7.2e-3,9.2e-3,11.2e-3};
	float trs[NPs] = {5e-3,10e-3,15e-3};
	float tips[NPs] = {6*PI/180,12*PI/180,80*PI/180};

	float TES[NACQS]{}; for (int i = 0; i < NACQS; i++) { TES[i] = tes[i%NTEs];}
	float TRS[NACQS]{}; for (int i = 0; i < NACQS; i++) { TRS[i] = trs[i/NTEs];}
	float TIPS[NACQS]{}; for (int i = 0; i < NACQS; i++) { TIPS[i] = tips[i/NTEs];}

	struct stopwatch sw;
	sw.click();
	float* pure_signal = noiseless_signal_vector(NACQS, TRS, TIPS, TES, 1.0, 312e-3, 822e-3, 50, 9050, 30, 0, 0, NUMTHREADS);
	sw.click();
	cout<<"Making orignal signal: "<<sw.check()<<endl;

	sw.click();
	float* pure_signal_mat = tile_vector(pure_signal,2*NACQS, NSIMS, NUMTHREADS);
	sw.click();
	cout<<"Repeating the Signal: "<<sw.check()<<endl;

	sw.click();
	float* simulated_signal_mat = add_noise2vector(pure_signal_mat, (2*NACQS*NSIMS), 1.0, NUMTHREADS);
	sw.click();
	cout<<"Creating Simulated Signal: "<<sw.check()<<endl;

	printmat(simulated_signal_mat, NSIMS, 2*NACQS);


	return 0;
}



