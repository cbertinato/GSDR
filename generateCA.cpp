// generateCA.cpp
#include "generateCA.h"

int op_mult (int i, int j) { return i*j; }

vector<int> generateCA (int sv) {
	// Code shift array
	// Shifts for the ground GPS transmitter are not included
	int g2s[] = {5,6,7,8,17,18,139,140,141,251, \
       				252,254,255,256,257,258,469,470,471,472,\
       				473,474,509,512,513,514,515,516,859,860,\
       				861,862,\
       				// end of shifts for GPS satellites
       				// Shifts for EGNOS and WAAS satellites (true_PRN = PRN + 87)
                 	145,175,52,21,237,235,886,657,\
       				634,762,355,1012,176,603,130,359,595,68,386};

	// G1 register
	vector<int> g1 (1023);
	deque<int> reg (10,-1);
	int readIn;

	// Generate all G1 signal chips based on the G1 feedback polynomial
	for (vector<int>::iterator i = g1.begin(); i <= g1.end(); ++i) {

		*i = reg.at(9);
		readIn = reg.at(2)*reg.at(9);
		reg.pop_front();
		reg.push_back(readIn);
	}

	// G2 register
	vector<int> g2 (1023);
	reg.assign(10,-1);

	// Generate all G2 signal chips based on the G2 feedback polynomial
	for (vector<int>::iterator i = g2.begin(); i <= g2.end(); ++i) {

		*i = reg.at(9);
		readIn = reg.at(1)*reg.at(2)*reg.at(5)*reg.at(7)*reg.at(8)*reg.at(9);
		reg.pop_front();
		reg.push_back(readIn);
	}

	// Shift G2
	int shift = g2s[sv];
	deque<int> g2Shift1 (g2.end() - shift + 1, g2.end());
	deque<int> g2Shift2 (g2.begin(), g2.end() - shift);

	// Binary add
	vector<int> CA (1023);

	transform( g1.begin(),g1.end(),g2.begin(),CA.begin(),op_mult );

	vector<int>::iterator i = CA.end();
	*i = op_mult( *g1.end(),*g2.end() );

	return CA;
}
