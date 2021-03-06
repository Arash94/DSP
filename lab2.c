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
#define wc PI/10

void convolve(float* h, float* Y, float* X);
void iir_notch(float y[], float x[], float c);
void sinc_filter(float* h, float c);
void add_signal(float* signal,float* noise,float* results);
void noise(float* noise, float c);
void echoed(float signal[],float output[]);
void normalizeArray(float* x);

int main(int argc, char *argv[])
{
	//Require 2 arguments: input file and output file
	//Require 2 arguments: input file and output file
	if(argc < 2)
	{
		printf("Not enough arguments \n");
		return -1;
	}


	SF_INFO sndInfo;
	SNDFILE *sndFile = sf_open(argv[1], SFM_READ, &sndInfo);

    
    SF_INFO sndInfoOut;
    sndInfoOut.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    sndInfoOut.channels = 1;
    sndInfoOut.samplerate = Fs;
    SNDFILE *sndFileOut = sf_open(argv[2], SFM_WRITE, &sndInfoOut);


	if (sndFile == NULL)
	{
		fprintf(stderr, "Error reading source file '%s': %s\n", argv[1], sf_strerror(sndFile));
		return 1;
	}

	float* hn = malloc(filter_size * sizeof(float));
    float* speech = malloc(theNumOfSamples * sizeof(float));
	float* echo = malloc(theNumOfSamples * sizeof(float));
    float* noise_signal = malloc(theNumOfSamples * sizeof(float));
    float* noisy_speech = malloc(theNumOfSamples * sizeof(float));
    float* filtered_sinc = malloc(theNumOfSamples * sizeof(float));
    float* filtered_iir = malloc(theNumOfSamples * sizeof(float));


	sf_readf_float(sndFile, speech, theNumOfSamples);
    
	//TASK 6: Add Echo

	echoed(speech, echo);
	sf_writef_float(sndFile, echo, theNumOfSamples);
	sf_write_sync(echo);
	free(echo);

	
	// //Taske 7: Filter FIR Notch
	// noise(noise_signal, wf);
	// add_signal(speech,noise_signal, noisy_speech);
	// sinc_filter(hn, wc);
	// convolve(hn, filtered_sinc,noisy_speech);
	// normalizeArray(filtered_sinc);
	// sf_writef_float(sndFile, filtered_sinc, theNumOfSamples);
	// free(filtered_sinc);
	//Task 8: Filter using IIR filter 
	// iir_notch(filtered_iir, noisy_speech, wc);
	// sf_writef_float(sndFile, filtered_iir, theNumOfSamples);
	// normalizeArray(filtered_iir);
	// free(filtered_iir);

	// sf_close(sndFile);
	// sf_write_sync(sndFile);
	// sf_close(sndFile);


	return 1;
}


void echoed(float signal[], float output[])
{
	float init= 0;

	printf("Creating the echoed version of the incoming signal\n");


	for(int i = 4000-1; i < theNumOfSamples; i++)
	{
		init = output[i-4000+1];
        output[i] = 0.7*output[i] + 0.3*init;
	}

}

void add_signal(float* signal,float* noise,float* results)
{   
    for (int i = 0; i < theNumOfSamples; i++)
    {
        *(results+i) = *(signal+i) + *(noise+i);
    }

}

void noise(float* noise, float c)
{
    for (int i = 0; i < theNumOfSamples; i++)
    {
       *(noise+i) = sin(2*PI*c*i);
    }   

}
void sinc_filter(float h[], float c)
{

    h[0] = c/PI;

    for (int i = 1; i < filter_size; i++)
    {
        *(h+i) = (c/PI*sin(c/PI*i))/(i*PI);
        if (i==7)
		{
			h[7] = h[7]-1.0;
		}
    }
   

}

void convolve(float* h, float* Y, float* X)
{
	for (int i = 7;i < theNumOfSamples-7;i++)
	{
       for (int j = 1; j < filter_size; j++)
	   {
           Y[i] = Y[i] + h[j]*X[i+j-7];
	   }
	}
}

void iir_notch(float y[], float x[], float c)
{
	//wc = pi/10
    float r1 = 1;
    float r2 = 0.95;
	float a[] = {1, -2*r1*cos(c), 1} ;
    float b[] = {1, -2*r2*cos(c), pow(r2,2)};

	y[0] = a[0] * x[0];
	y[1] = a[0] * x[1] + a[1] * x[0] - b[1] * y[0];
	// y[2] = a[0] * x[3] + a[1] * x[2] + a[2] * y[1] - b[1] * y[1] - b[2] * y[0];

	for(int i = 2; i< theNumOfSamples; i++)
	{
		y[i] = a[0] * x[i] + a[1] * x[i-1] + a[2] * x[i-2] - b[1] * y[i-1] - b[2] * y[i-2];
	}
}

void normalizeArray(float x[])
{
	float maxAbsVal = 0.0;
	for(int i = 0; i < theNumOfSamples; i++)
	{
		if(fabs(x[i]) > maxAbsVal)
		{
			maxAbsVal = fabs(x[i]);
		}
	}
	for(int i = 0; i<theNumOfSamples; i++)
	{
		x[i] = x[i]/maxAbsVal;
	}
}
