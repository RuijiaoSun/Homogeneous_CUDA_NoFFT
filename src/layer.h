#ifndef LAYER_H
#define LAYER_H
#include <vector>
#include <complex>
using namespace std;

void crossProduct(vector<complex<double>>* a, vector<complex<double>>* b, vector<complex<double>>* c);
complex<double> Norm(vector<complex<double>>* E);
void Normalize(vector<complex<double>>* E);
vector<vector<double> > transpose(vector<vector<double> >* matrix);
void orthogonalize(vector<complex<double>>* E_0, vector<complex<double>>* E0, vector<complex<double>>* d);

class layersample {

public:
	double k;						//wavenumber.
	vector<complex<double>> n;		//refractive index.	
	vector<double> z;				//z postions.
	vector<complex<double>> s;		//propagation direction. Keep it for output.
	vector<complex<double>> sz;		//propagation direction. Keep it for output.
	vector<complex<double>> d;		//direction of propagation of the plane wave.
	vector<complex<double>> Pt;		//transimission.
	vector<complex<double>> Pr;		//reflection.
	//Calulate the index of the field component associated with a layer.
	// l is the layer index. c is the component(x, y, z). d is the direction(0for transmission, 1 for reflection).
	size_t ii(size_t l, size_t c, size_t d);

	//Generate the linear system corresponding to this layered sample and plane wave.
	void generate_linsys(size_t LAYERS, vector<complex<double>>& M, vector<complex<double>>& b, vector<complex<double>>& E, vector<complex<double>>* P, bool CPU_op);

	//Build matrix and get E.
	void solve(vector<complex<double>>* E, bool CPU_op);
};

#endif