// acquire.cpp

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <fftw3.h>		// Fast Fourier Transform

using namespace std;

# define TRUE 1
# define FALSE 0

const double PI = 3.141592;

// Declarations
struct settings_t {
	double SF;							// Sampling frequency
	double IF;							// Intermediate frequency
	double CF;							// Code frequency
	unsigned int codeLength;			// Code length in chips
	double acqSearchBand;				// Frequency search bandwidth
	unsigned long int samplesPerCode;
	bool GPS;
	bool GLONASS;
	unsigned int acqFreqStep;
	unsigned int acqBandwidth;
} settings;

typedef struct acqData_t {
	bool acquired;
	double carrierFreq;
	unsigned int codeOffset;
	signed int doppler;
	float acqMetric;
} acqData;

// Function declarations
double mean (vector<double>);
vector<int> generateCA (int);
vector<double>::iterator max (vector<double>);
fftw_complex * conj(fftw_complex*, unsigned int);

int main(int argc, char* argv[]) {

	if (argc < 2) {
		cerr << "Usage: " << argv[0] << "<data file path>" << endl;
		return 1;
	}

	ifstream data;

	char *inputFile = argv[1];

	settings.SF = 38.192e6; // Hz
	settings.IF = 9.548e6;  // Hz
	settings.CF = 1.023e6;  // Hz
	settings.codeLength = 1023; // chips
	settings.samplesPerCode = round(settings.SF / (settings.CF / settings.codeLength));
	settings.GPS = TRUE;
	settings.GLONASS = FALSE;
	settings.acqFreqStep = 500;		// Hz
	settings.acqBandwidth = 20e3;	// Hz

	vector<int> satList (32);

	if (settings.GPS && settings.GLONASS)
		satList.resize(55);
	else if (!settings.GPS && !settings.GLONASS) {
		cout << "No satellites!" << endl;
		return -1;
	}

	int i = 1;
	for ( vector<int>::iterator t = satList.begin(); t != satList.end(); ++t )
		*t = i++;

	cout << "-------------------------------" << endl;
	cout << "Sample Frequency: " << settings.SF << " Hz" << endl;
	cout << "Intermediate Frequency: " << settings.IF << " Hz" << endl;
	cout << "Code Frequency Basis: " << settings.CF << " Hz" << endl;
	cout << "Code Length: " << settings.codeLength << endl;
	cout << "Samples per Code: " << settings.samplesPerCode << endl;
	cout << "-------------------------------" << endl;

	data.open(inputFile,ios::binary);

	if (!data) {
		cout << "Error opening input data file." << endl;
		return -1;
	}

	vector<char> signal (settings.samplesPerCode*11);

	// Read 11 ms of data
	try {
		data.read((char*)&signal[0],settings.samplesPerCode*sizeof(char)*11);
	} catch (int e) {
		cout << "An exception occurred while reading data. (" << e << ")" << endl;
	}

	cout << "Read " << data.gcount() << " bytes." << endl;

	// Create two 1 msec signal vectors
	vector<double> signal1 (settings.samplesPerCode);
	vector<double> signal2 (settings.samplesPerCode);

	for (vector<char>::size_type j = 0; j < settings.samplesPerCode; j++)
	{
		signal1[j] = (double)signal[j];
		signal2[j] = (double)signal[j + settings.samplesPerCode];

		//cout << "Signal: " << signal[j] << ", Signal1: " << signal1[j] << endl;
	}

	// Create zero DC signal vector
	vector<char> signal0DC (settings.samplesPerCode);

	// Compute the number of frequency bins for the given search band; steps of 500 Hz
	// (IF +/- acqBandwidth)/500 Hz + 1
	// Number of frequencies to iterate through
  unsigned int numFrequencies = round(settings.acqBandwidth / settings.acqFreqStep) + 1;
	cout << "Number of frequency bins: " << numFrequencies << endl;

	cout << "Satellite acquisition list: ";
	for ( vector<int>::iterator t = satList.begin(); t!=satList.end(); ++t) {
		cout << *t << " ";
	}
	cout << endl;

	vector<acqData> satAcqData (satList.size());
	vector<acqData>::iterator j = satAcqData.begin();

  // Initialize variables for fft
  fftw_plan p1, p2;
  fftw_complex *sampledPRN, *sampledPRNfft, *X1, *X2, *X1fft, *X2fft, *temp, *tempIFFT;

  sampledPRN = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * settings.samplesPerCode);
  sampledPRNfft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * settings.samplesPerCode);
	X1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * settings.samplesPerCode);
	X2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * settings.samplesPerCode);
	X1fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * settings.samplesPerCode);
	X2fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * settings.samplesPerCode);
	temp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * settings.samplesPerCode);
	tempIFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * settings.samplesPerCode);

	// Iterate through satellites
	for (vector<int>::iterator i = satList.begin(); i != satList.end(); ++i) {

		cout << "Searching for SV " << *i << "..." << endl;

		// Generate PRN
		vector<int> PRN = generateCA(*i);
		vector<int>::iterator j = PRN.begin();

		p1 = fftw_plan_dft_1d(settings.samplesPerCode, sampledPRN, sampledPRNfft, FFTW_FORWARD, FFTW_ESTIMATE);

		// Sample PRN
		unsigned int n = 1;
		for (unsigned int k = 0; k < settings.samplesPerCode; k++) {
        	sampledPRN[k][0] = *(j++);
        	sampledPRN[k][1] = 0;
        	n = floor( (k+1)*settings.CF / settings.SF ) + 1;
		}

		// Compute fft of sampled PRN signal
		fftw_execute(p1);
		//sampledPRNfft = conj(sampledPRNfft, settings.samplesPerCode);

		// Loop through frequency band
		double carrierFreq;
		vector<double> Y1 (settings.samplesPerCode);
		vector<double> Y2 (settings.samplesPerCode);
		vector<double> maxPower (numFrequencies);
		vector<double>::iterator mp = maxPower.begin();

		for (unsigned int k = 0; k < numFrequencies; k++) {

			// Generate carrier signal and wipe off
			carrierFreq = settings.IF - 10e3 + settings.acqFreqStep*k;

			p1 = fftw_plan_dft_1d(settings.samplesPerCode, X1, X1fft, FFTW_FORWARD, FFTW_ESTIMATE);
			p2 = fftw_plan_dft_1d(settings.samplesPerCode, X2, X2fft, FFTW_FORWARD, FFTW_ESTIMATE);

			for (unsigned int n = 0; n < settings.samplesPerCode; n++) {
				X1[n][0] = signal1[n]*cos(2*PI*carrierFreq*1/settings.SF*n);
				X1[n][1] = signal1[n]*sin(2*PI*carrierFreq*1/settings.SF*n);

				//cout << "Signal1[" << n << "] = " << signal1[n] << endl;
				//cout << "k=" << k << "; X1[" << n << "] = " << X1[n][0] << " + i" << X1[n][1] << endl;

				X2[n][0] = signal2[n]*cos(2*PI*carrierFreq*1/settings.SF*n);
				X2[n][1] = signal2[n]*sin(2*PI*carrierFreq*1/settings.SF*n);
			} // carrier wipe-off loop

			// FFT
			fftw_execute(p1); // ok?
			fftw_execute(p2); // ok?

			//cout << "k=" << k << "; X1fft[1] = " << X1fft[1][0] << " + i" << X1fft[1][1] << endl;

			p1 = fftw_plan_dft_1d(settings.samplesPerCode, temp, tempIFFT, FFTW_BACKWARD, FFTW_ESTIMATE);

			// Power
			for (unsigned int n = 0; n < settings.samplesPerCode; n++) {
				temp[n][0] = sampledPRNfft[n][0]*X1fft[n][0] - sampledPRNfft[n][1]*X1fft[n][1];
				temp[n][1] = sampledPRNfft[n][0]*X1fft[n][1] + sampledPRNfft[n][1]*X1fft[n][0];
			}

			// inverse FFT
			fftw_execute(p1);

			p1 = fftw_plan_dft_1d(settings.samplesPerCode, temp, tempIFFT, FFTW_BACKWARD, FFTW_ESTIMATE);

			for (unsigned int n = 0; n < settings.samplesPerCode; n++) {
				Y1[n] = tempIFFT[n][0]*tempIFFT[n][0] + tempIFFT[n][1]*tempIFFT[n][1];

				//cout << "tempIFFT = " << tempIFFT[n][0] << " + i" << tempIFFT[n][1] << endl;
				//cout << "Y1: " << Y1[n] << endl;

				temp[n][0] = sampledPRNfft[n][0]*X2fft[n][0] - sampledPRNfft[n][1]*X2fft[n][1];
				temp[n][1] = sampledPRNfft[n][0]*X2fft[n][1] + sampledPRNfft[n][1]*X2fft[n][0];

			}

			// inverse FFT
			fftw_execute(p1);

			for (unsigned int n = 0; n < settings.samplesPerCode; n++)
				Y2[n] = tempIFFT[n][0]*tempIFFT[n][0] + tempIFFT[n][1]*tempIFFT[n][1];

			// Check which msec had the greatest power and save it for each frequency
			vector<double>::iterator it1 = max(Y1);
			vector<double>::iterator it2 = max(Y2);

			//cout << "Max Y1: " << *it1 << endl;
			//cout << "Max Y2: " << *it2 << endl;

			if (*it1 > *it2)
				*(mp++) = *it1;
			else
				*(mp++) = *it1;

		} // frequency loop

		// Search for correlation peak
		vector<double>::iterator corr = max(maxPower);
		double corrFreq = distance(maxPower.begin(),corr);
		cout << "Max Power: " << *corr << endl;
		cout << "Frequency Index: " << corrFreq << endl;

	}

	// Clean up
	data.close();
	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);
	fftw_free(sampledPRN);
	fftw_free(sampledPRNfft);

	fftw_free(temp);
	fftw_free(tempIFFT);

	fftw_free(X1);
	fftw_free(X2);
	fftw_free(X1fft);
	fftw_free(X2fft);

	return 0;
}
