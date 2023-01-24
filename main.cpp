#include <iostream>
#include <complex>
#define _USE_MATH_DEFINES
#include "math.h"
#include <time.h>
#include <memory>

using namespace std;

const int MAX_ORDER = 13;
const bool PRINT_COEFS = false;
const std::complex<double> constI(0, 1);
const double constE = exp(1.0);


template <typename A>
complex<A>* dft(A arr[], const int arrSize){
    complex<A>* result = new complex<A>[arrSize];
    for(int i = 0; i < arrSize; i++){
        double rel = 0;
        double com = 0;
        for(int j = 0; j < arrSize; j++){
            double alpha = 2 * M_PI * i * j / arrSize;
            rel += arr[j] * cos(alpha);
            com += arr[j] * sin(alpha);
        }
        complex<A> temp(rel, com);
        result[i] = temp;
    }
    return result;
}

template <typename B>
complex<B>* fft(B arr[], const int arrSize){
    complex<B>* store = new complex<B>[arrSize];

    for(int i = 0; i < arrSize; i++){
        double alpha = -2 * M_PI * i / arrSize;
        double rel = cos(alpha), com = sin(alpha);
        complex<B> temp(rel, com);
        store[i] = temp;
    }

    complex<B>* x0 = new complex<B>[arrSize/2];
    complex<B>* x1 = new complex<B>[arrSize/2];

    for(int i = 0; i < arrSize/2; i++){
        x0[i] = arr[i*2];
        x1[i] = arr[i*2 + 1];
    }

    complex<B>* y0 = fft(x0, arrSize/2);
    complex<B>* y1 = fft(x1, arrSize/2);

    complex<B>* result = new complex<B>[arrSize];

    for(int i = 0; i < arrSize/2; i++){
        result[i] = y0[i] + store[i] * y1[i];
        result[i + arrSize / 2] = y0[i] - store[i] * y1[i];
    }
    return result;
}

int main()
{
//    for(int o = 1; o <= MAX_ORDER; o++){
//        const int N = 1 << o;
//        printf("N: %i\n", N);
//
//        double* f = new double[N];
//        for(int n = 0; n < N; n++){
//            f[n] = n / (double)N;
//        }
//
//        clock_t t1 = clock();
//        complex<double>* cDFT = dft(f, N);
//        clock_t t2 = clock();
//        double dft_time = (t2 - t1) / (double)CLOCKS_PER_SEC * 1000.0;
//        printf("DFT time [ms]: %f\n", dft_time);
//
//        t1 = clock();
//        complex<double>* cFFT = fft(f, N);
//        t2 = clock();
//        double fft_time = (t2 - t1) / (double)CLOCKS_PER_SEC * 1000.0;
//        printf("FFT time [ms]: %f\n", dft_time);
//
//        printf("mean error: %f\n", err(cDFT, cFFT, N));
//
//        if(PRINT_COEFS){
//            for(int k = 0; k < N; k++){
//                  cout << "x" << k << ": DFT " << cDFT[k] << " FFT " << cFFT[i] << endl;
//            }
//        }
//        delete[] f;
//        delete[] cDFT;
//        delete[] cFFT;
//    }
//
//    return 0;

    const int N = 1 << 2;

    double* f = new double[N];
    for(int n = 1; n < N; n++){
        f[n] = n / (double)N;
    }
    complex<double>* cFFT = fft(f, N);

    for(int i = 1; i < N; i++){
        cout << f[i] << ": ( " << cFFT[i].real() << ", " << cFFT[i].imag() << "i )" << endl;
    }

    delete[] cFFT;

    return 0;
}
