#include <stdlib.h>
#include <stdio.h>
#include <float.h>
//#include "wave.h"
#include <sndfile.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265
#define Fs 48000
#define theNumOfSamples 960000
#define filter_size 15
#define wf 2400

void convolve(float* Y, float *X, float* h);
void iir(float* y, float* x);
void sinc_filter(float* h, float wc);
void add_signal(float* signal,float* noise,float* results);
void noise(float* noise, float wf);


int main(int argc, char *argv[])
{
	//Require 2 arguments: input file and output file
	if(argc < 2)
	{
		printf("Not enough arguments \n");
		return -1;
	}

	SF_INFO sndInfoOut;
	sndInfoOut.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
	sndInfoOut.channels = 1;
	sndInfoOut.samplerate = Fs;
	SNDFILE *sndFileOut = sf_open(argv[1], SFM_WRITE, &sndInfoOut);
    
	 float* h = malloc(filter_size * sizeof(float));
     float* speech = malloc(theNumOfSamples * sizeof(float));
     float* noise = malloc(theNumOfSamples * sizeof(float));
     float* noisy_speech = malloc(theNumOfSamples * sizeof(float));
     float* filtered_since = malloc(theNumOfSamples * sizeof(float));
     float* filtered_iir = malloc(theNumOfSamples * sizeof(float));




	return 1;
}

void add_signal(float* signal,float* noise,float* results)
{   
    for (float i = 0; i < theNumOfSamples; i++)
    {
        results[i] = signal[i] + noise[i];
    }

}

void noise(float* noise, float wf)
{
    for (float i = 0; i < theNumOfSamples; i++)
    {
       noise[i] = sin(2*PI*wf*i);
    }   

}
void sinc_filter(float* h, float wc)
{

    h[0] = wc*/PI;

    for (float i = 1; i < filter_size; i++)
    {
        h[i] = (wc*/PI*sin(wc*/PI*i))/(i*PI);
        if (i==7){h[7] = h[7]-1.0};
    }
   

}

void convolve(float* h, float* Y, float* X)
{
	for (i = 7;i < theNumOfSamples-7;i++){
       for (j = 1; j < filter_size; j++)
	   {
           Y(i) = Y(i) + h(j)*X(i+j-7);
	   }
	}

}

void iir_notch(float* y, float* x)
{
	//wc = pi/10
    float r1 = 1;
    float r2 = 0.95;
	float a = {1, -2*r1*cos(wc), 1} ;
    float b = {1, -2*r2*cos(wc), r2**2};

	y[0] = a[0] * x[0];
	y[1] = a[0] * x[1] + a[1] * x[0] - b[1] * y[0];
	y[2] = a[0] * x[3] + a[1] * x[2] + a[2] * y[1] - b[1] * y[1] - b[2] * y[0];

	for(int i = 5; i< theNumOfSamples; i++){
		y[i] = a[0] * x[i] + a[1] * x[i-1] + a[2] * x[i-2] - b[1] * y[i-1] - b[2] * y[i-2];
	}
}