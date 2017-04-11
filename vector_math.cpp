#include "vector_math.h"

double mean (vector<double> a) {
	double b = 0;
	vector<double>::iterator i;

	for (i = a.begin(); i != a.end(); ++i)
	{
		b += *i;
	}

	b = b/a.size();

	return b;
}

vector<double>::iterator max (vector<double> a) {
	vector<double>::iterator max = a.begin();
	vector<double>::iterator i;
	//int length = a.size();
	//double max = a[0];

	for(i = a.begin(); i < a.end(); ++i)
	{
		if (*i > *max)
			max = i;
	}

	return max;

}

fftw_complex * conj(fftw_complex *c, unsigned int N) {

	for (unsigned int i=0; i < N; i++) {
		c[i][1] = -c[i][1];
	}

	return c;
}

fftw_complex * conj(fftw_complex *c) {

	c[0][1] = -c[0][1];
	return c;
}
