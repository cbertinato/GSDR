// fft.c
#include <complex>
#include <fftw3.h>

complex<double> * fft(complex<double> *a, unsigned int N) {
	fftw_complex *in, *b;
	b = reinterpret_cast<fftw_complex*>(a);
	
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	memcpy( &in, &b, sizeof(fftw_complex) * N );
}

complex<double> * ifft(*complex<double> a, unsigned int N) {
	fftw_complex *in, *b;
	b = reinterpret_cast<fftw_complex*>(a);
	
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	memcpy( &in, &b, sizeof(fftw_complex) * N );
}
