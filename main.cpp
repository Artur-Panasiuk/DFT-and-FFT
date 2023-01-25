#include <iostream>
#include <complex>
#define _USE_MATH_DEFINES
#include "math.h"
#include <time.h>
#include <memory>

using namespace std;

const int MAX_ORDER = 5;
const bool PRINT_COEFS = true;
const std::complex<double> constI(0, 1);
const double constE = exp(1.0);


complex<double>* dft(double arr[], const int arrSize) {
    complex<double>* result = new complex<double>[arrSize];
    for (int i = 0; i < arrSize; i++) {
        double rel = 0;
        double com = 0;
        for (int j = 0; j < arrSize; j++) {
            double alpha = 2 * M_PI * i * j / arrSize;
            rel += arr[j] * cos(alpha);
            com -= arr[j] * sin(alpha);
        }
        complex<double> temp(rel, com);
        result[i] = temp;
    }
    return result;
}


complex<double>* fft(double arr[], const int arrSize) {

    complex<double>* result = new complex<double>[arrSize];

    if (arrSize == 1) {

        for (int i = 0; i < arrSize; i++) {
            result[i] = arr[i];
        }

        return result;
    }

    complex<double>* store = new complex<double>[arrSize];

    for (int i = 0; i < arrSize; i++) {
        double alpha = -2 * M_PI * i / arrSize;
        double rel = cos(alpha), com = sin(alpha);
        complex<double> temp(rel, com);
        store[i] = temp;
    }

    double* x0 = new double[arrSize / 2];
    double* x1 = new double[arrSize / 2];

    for (int i = 0; i < arrSize / 2; i++) {
        x0[i] = arr[i * 2];
        x1[i] = arr[i * 2 + 1];
    }

    complex<double>* y0 = fft(x0, arrSize / 2);
    complex<double>* y1 = fft(x1, arrSize / 2);

    for (int i = 0; i < arrSize / 2; i++) {
        result[i] = y0[i] + store[i] * y1[i];
        result[i + arrSize / 2] = y0[i] - store[i] * y1[i];
    }
    delete[] store;
    delete[] x0;
    delete[] x1;

    return result;
}

double err(complex<double> dftArr[], complex<double> fftArr[], int N) {
    double result = 0;
    for (int i = 0; i < N; i++) {
        complex<double> tempComplex = dftArr[i] - fftArr[i];
        result += abs(tempComplex);
    }
    return result / N;
}

int main()
{
        for(int o = 1; o <= MAX_ORDER; o++){
            const int n = 1 << o;
            printf("n: %i\n", n);
    
            double* f = new double[n];
            for(int i = 0; i < n; i++){
                f[i] = i / (double)n;
            }
    
            clock_t t1 = clock();
            complex<double>* cdft = dft(f, n);
            clock_t t2 = clock();
            double dft_time = (t2 - t1) / (double)CLOCKS_PER_SEC * 1000.0;
            printf("dft time [ms]: %f\n", dft_time);
    
            t1 = clock();
            complex<double>* cfft = fft(f, n);
            t2 = clock();
            double fft_time = (t2 - t1) / (double)CLOCKS_PER_SEC * 1000.0;
            printf("fft time [ms]: %f\n", fft_time);
    
            printf("mean error: %lg\n", err(cdft, cfft, n));
    
            if(PRINT_COEFS){
                for(int k = 0; k < n; k++){
                      cout << "x" << k << ": dft " << cdft[k] << " fft " << cfft[k] << endl;
                }
            }
            delete[] f;
            delete[] cdft;
            delete[] cfft;
        }
    
        return 0;
}
