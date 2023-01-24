#include <iostream>
#include <complex>
#include "math.h"
#include <time.h>

using namespace std;

const int MAX_ORDER = 13;
const bool PRINT_COEFS = false;


template <typename A>
complex<A>* dft(A arr[], const int arrSize){
    double r[arrSize];
    double j[arrSize];
    for(int i = 0; i < arrSize; i++){

    }
}

template <typename B>
complex<B>* fft(B arr[], int arrSize){

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
//                //wypis na ekran wspolczynnikow obu transformat ( czesci rzeczywiste i urojone )
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
    for(int n = 0; n < N; n++){
        f[n] = n / (double)N;
        cout << f[n]<<endl;
    }
    complex<double>* cDFT = dft(f, N);

    delete[] cDFT;

    return 0;
}
