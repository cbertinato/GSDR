// vector_math.h

#ifndef VECTOR_MATH_H
#define VECTOR_MATH_H

#include <vector>
#include <fftw3.h>

using namespace std;

double mean (vector<double>);
vector<double>::iterator max (vector<double>);
fftw_complex * conj(fftw_complex *, unsigned int);
fftw_complex * conj(fftw_complex *);

#endif
