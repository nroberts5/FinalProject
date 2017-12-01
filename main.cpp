#include <iostream>
#include <cstring>
#include <cmath>
#include "H5Cpp.h"

using namespace std;

#define FIELD_STRENGTH 3.0

float* signal_equation(const float TR, const float alpha, const float TE, const float beta, const float T1__F, const float T1__W, const float rho__F, const float rho__W, const float R2s, const float phi, const float psi){
	float ppm[6] = {0.6, -0.5, -1.95, -2.60, -3.40, -3.80};
	float fp[6]; memcpy(fp,ppm,6*sizeof(float)); for(int i = 0; i < 6; i++){fp[i] *= (FIELD_STRENGTH*42.577);}
    float a[6] = {4.7/ 100, 3.9/ 100, 0.6/ 100, 12.0/ 100, 70.0/ 100, 8.8/ 100};
	
	float sR = ((rho__W * sin(beta * alpha) * (0.1e1 - exp(-TR / T1__W)) / (0.1e1 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (0.1e1 - exp(-TR / T1__F)) * (a[0] * cos(0.2e1 * 0.3141592654e1 * fp[0] * TE) + a[1] * cos(0.2e1 * 0.3141592654e1 * fp[1] * TE) + a[2] * cos(0.2e1 * 0.3141592654e1 * fp[2] * TE) + a[3] * cos(0.2e1 * 0.3141592654e1 * fp[3] * TE) + a[4] * cos(0.2e1 * 0.3141592654e1 * fp[4] * TE) + a[5] * cos(0.2e1 * 0.3141592654e1 * fp[5] * TE)) / (0.1e1 - cos(beta * alpha) * exp(-TR / T1__F))) * cos(TE * psi + phi) - rho__F * sin(beta * alpha) * (0.1e1 - exp(-TR / T1__F)) * (a[0] * sin(0.2e1 * 0.3141592654e1 * fp[0] * TE) + a[1] * sin(0.2e1 * 0.3141592654e1 * fp[1] * TE) + a[2] * sin(0.2e1 * 0.3141592654e1 * fp[2] * TE) + a[3] * sin(0.2e1 * 0.3141592654e1 * fp[3] * TE) + a[4] * sin(0.2e1 * 0.3141592654e1 * fp[4] * TE) + a[5] * sin(0.2e1 * 0.3141592654e1 * fp[5] * TE)) * sin(TE * psi + phi) / (0.1e1 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
	float sI = ((rho__W * sin(beta * alpha) * (0.1e1 - exp(-TR / T1__W)) / (0.1e1 - cos(beta * alpha) * exp(-TR / T1__W)) + rho__F * sin(beta * alpha) * (0.1e1 - exp(-TR / T1__F)) * (a[0] * cos(0.2e1 * 0.3141592654e1 * fp[0] * TE) + a[1] * cos(0.2e1 * 0.3141592654e1 * fp[1] * TE) + a[2] * cos(0.2e1 * 0.3141592654e1 * fp[2] * TE) + a[3] * cos(0.2e1 * 0.3141592654e1 * fp[3] * TE) + a[4] * cos(0.2e1 * 0.3141592654e1 * fp[4] * TE) + a[5] * cos(0.2e1 * 0.3141592654e1 * fp[5] * TE)) / (0.1e1 - cos(beta * alpha) * exp(-TR / T1__F))) * sin(TE * psi + phi) + rho__F * sin(beta * alpha) * (0.1e1 - exp(-TR / T1__F)) * (a[0] * sin(0.2e1 * 0.3141592654e1 * fp[0] * TE) + a[1] * sin(0.2e1 * 0.3141592654e1 * fp[1] * TE) + a[2] * sin(0.2e1 * 0.3141592654e1 * fp[2] * TE) + a[3] * sin(0.2e1 * 0.3141592654e1 * fp[3] * TE) + a[4] * sin(0.2e1 * 0.3141592654e1 * fp[4] * TE) + a[5] * sin(0.2e1 * 0.3141592654e1 * fp[5] * TE)) * cos(TE * psi + phi) / (0.1e1 - cos(beta * alpha) * exp(-TR / T1__F))) * exp(-R2s * TE);
	
	float* sig = new float[2];
	sig[0] = sR;
	sig[1] = sI;
	return sig;
}


int main(int argc, char const *argv[])
{
  /* code */
  float* test = signal_equation(20e-3,60*3.14/180,1.2e-3,1.0,312e-3,822e-3,50,9050,30,0,0);
  cout<<test[0]<<"+j"<<test[1]<<endl;



  return 0;
}