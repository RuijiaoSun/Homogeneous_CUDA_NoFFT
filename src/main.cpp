//Update the data output format in example_layer().
#include "layer.h"
#include <windows.h>

#include <fstream>
#include <iostream>
#include "stim/parser/arguments.h"

using namespace std;
#define PI 3.14159265358979323846264338328

bool ASCII_output = false;
bool CPU_op = false;

void advertise() {
	std::cout << std::endl << std::endl;
	std::cout << "=========================================================================" << std::endl;
	//std::cout << "Thank you for using the NetMets network comparison tool!" << std::endl;
	std::cout << "Scalable Tissue Imaging and Modeling (STIM) Lab, University of Houston" << std::endl;
	std::cout << "=========================================================================" << std::endl << std::endl;
}

//Calculate Ez from vector d and Ex Ey.
void E_Cal(vector<complex<double>>* E, vector<complex<double>>* d) {
	if (d->size() == 2) {
		complex<double> dz = sqrt(1.0 + 0i - pow((*d)[0], 2) - pow((*d)[1], 2));
		d->push_back(dz);
	}
	E->push_back(-((*E)[0] * (*d)[0] + (*E)[1] * (*d)[1]) / (*d)[2]);
}
void output_binary(std::string filename, layersample& layers) {
	size_t L = layers.n.size();											// get the number of layers
	std::ofstream outFile;
	outFile.open(filename, std::ofstream::binary);						// open the output file for binary writing
	if (outFile) {
		outFile.write((char*)&layers.k, sizeof(double));
		outFile.write((char*)&layers.d[0], sizeof(double));
		outFile.write((char*)&layers.d[1], sizeof(double));
		outFile.write((char*)&layers.n[0], sizeof(double) * 2);

		for (size_t i = 0; i < L; i++) {
			outFile.write((char*)&layers.z[i], sizeof(double));
			outFile.write((char*)&layers.sz[i], 2 * sizeof(double));
			outFile.write((char*)&layers.Pt[i], 2 * sizeof(double));
			outFile.write((char*)&layers.Pt[i + L], 2 * sizeof(double));
			outFile.write((char*)&layers.Pt[i + 2 * L], 2 * sizeof(double));
			outFile.write((char*)&layers.Pr[i], 2 * sizeof(double));
			outFile.write((char*)&layers.Pr[i + L], 2 * sizeof(double));
			outFile.write((char*)&layers.Pr[i + 2 * L], 2 * sizeof(double));
		}
		outFile.close();
	}
	else {
		std::cout << "ERROR opening output file for binary writing: " << filename << std::endl;
	}
}

void output_txt(std::string filename, layersample& layers) {
	size_t L = layers.n.size();											// get the number of layers
	std::ofstream outFile;
	outFile.open(filename);												// open the output file for text writing
	if (!outFile) {
		std::cout << "ERROR: Could not open file " << filename << std::endl;
		exit(1);
	}
	int width = 15;

	outFile << "--------------------------------" << endl;
	outFile << "The wavenumber at free space is : " << layers.k << endl;
	for (size_t i = 0; i < L; i++) {
		if (i == 0) {
			outFile << "--------------------------------" << endl;
			outFile << "LAYER " << i << " (z = " << layers.z[i] << ")" << endl;
			outFile << "refractive index: " << layers.n[i].real() << " + i " << layers.n[i].imag() << endl;
			outFile << "----------------------" << endl;
			outFile << "sx = " << setw(width) << layers.s[0].real() << " + i " << layers.s[0].imag() << endl;
			outFile << "sy = " << setw(width) << layers.s[1].real() << " + i " << layers.s[1].imag() << endl;
			outFile << "sz = " << setw(width) << layers.sz[i].real() << " + i " << layers.sz[i].imag() << endl;
			outFile << "----->>>>>" << endl;
			outFile << " X = " << setw(width) << layers.Pt[i].real() << " + i " << layers.Pt[i].imag() << endl;
			outFile << " Y = " << setw(width) << layers.Pt[i + L].real() << " + i " << layers.Pt[i + L].imag() << endl;
			outFile << " Z = " << setw(width) << layers.Pt[i + 2 * L].real() << " + i " << layers.Pt[i + 2 * L].imag() << endl;
			outFile << "<<<<<-----" << endl;
			outFile << " X = " << setw(width) << layers.Pr[i].real() << " + i " << layers.Pr[i].imag() << endl;
			outFile << " Y = " << setw(width) << layers.Pr[i + L].real() << " + i " << layers.Pr[i + L].imag() << endl;
			outFile << " Z = " << setw(width) << layers.Pr[i + 2 * L].real() << " + i " << layers.Pr[i + 2 * L].imag() << endl;
		}
		else {
			outFile << "----------------------" << endl;
			outFile << "LAYER " << i << " (z = " << layers.z[i] << ")" << endl;
			outFile << "refractive index: " << layers.n[i].real() << " + i " << layers.n[i].imag() << endl;
			outFile << "----------------------" << endl;
			outFile << "sx = " << setw(width) << layers.s[0].real() << " + i " << layers.s[0].imag() << endl;
			outFile << "sy = " << setw(width) << layers.s[1].real() << " + i " << layers.s[1].imag() << endl;
			outFile << "sz = " << setw(width) << layers.sz[i].real() << " + i " << layers.sz[i].imag() << endl;
			outFile << "----->>>>>" << endl;
			outFile << " X = " << setw(width) << layers.Pt[i].real() << " + i " << layers.Pt[i].imag() << endl;
			outFile << " Y = " << setw(width) << layers.Pt[i + L].real() << " + i " << layers.Pt[i + L].imag() << endl;
			outFile << " Z = " << setw(width) << layers.Pt[i + 2 * L].real() << " + i " << layers.Pt[i + 2 * L].imag() << endl;

			outFile << "<<<<<-----" << endl;
			outFile << " X = " << setw(width) << layers.Pr[i].real() << " + i " << layers.Pr[i].imag() << endl;
			outFile << " Y = " << setw(width) << layers.Pr[i + L].real() << " + i " << layers.Pr[i + L].imag() << endl;
			outFile << " Z = " << setw(width) << layers.Pr[i + 2 * L].real() << " + i " << layers.Pr[i + 2 * L].imag() << endl;
		}

	}
	outFile.close();
}

void calculate_layer(std::string outName,
					vector<complex<double>>* ns,		// the refractive index.
					vector<double>* depths,				// z-direction position.
					vector<complex<double>>* E0,		// the initialized E0.
					vector<complex<double>>* d0,			// direction of propagation of the plane wave.
					double k0) {						// the wavenumber at free space.
	std::string outName_ext(outName, size(outName)-4, size(outName)-1);			// extract the extension of the output file				
	E_Cal(E0, d0);					// make sure that both vectors are orthogonal.

	//Creat a new layersample and initialize.
	layersample Layer1;									// create a layered sample
	Layer1.n = *ns;										// set a pointer to the refractive indices
	Layer1.z = *depths;									// set a pointer to the layer depths
	Layer1.d = *d0;
	Layer1.k = k0;

	LARGE_INTEGER t1, t2, tc;							// Timing.
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
	Layer1.solve(E0, CPU_op);			// Solve for the substrate field in GPU.
	QueryPerformanceCounter(&t2);
	std::cout << "time for 'solving linear functions':" << (t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart << "ms."<< std::endl;
	//output(outName, Layer1);
	if (ASCII_output)
		output_txt(outName, Layer1);		
	else
		output_binary(outName, Layer1);
}

//Main function for example_layer.
int main(int argc, char* argv[]) {
	stim::arglist args;

	//Basic argument lists.
	args.add("help", "prints this help");
	args.add("s", "propagation direction vector (x, y)", "0.5 0.0", "[-1.0, 1.0]");
	args.add("l", "wavelength", "5.0", "in arbitrary units (ex. um)");
	args.add("Ex", "complex amplitude (x direction)", "0.5 0.0");
	args.add("Ey", "complex amplitude (y direction)", "0.5 0.0");
	args.add("z", "layer positions");
	args.add("n", "layer optical path length (real refractive index)", "1.0 1.4 1.4 1.0");
	args.add("kappa", "layer absorbance (imaginary refractive index)");
	args.add("ascii", "output as an ASCII file");
	args.add("CPU", "execute the program in CPU");
	args.parse(argc, argv);

	if (args["help"].is_set()) {				// test for help
		advertise();							// output the advertisement
		std::cout << args.str();				// output arguments
		exit(1);								// exit
	}

	ASCII_output = args["ascii"];			// if the ascii flag is set, output an ascii file
	CPU_op = args["CPU"];					// if the CPU_op is set, run the program in CPU.
	std::string outName;
	if (args.nargs() == 1) {
		outName = args.arg(0);
	}
	else if (args.nargs() == 0) {
		if (ASCII_output)
			outName = "output.txt";
		else
			outName = "output.lyr";
	}
	else {
		std::cout << "ERROR: Too many arguments." << std::endl;
		exit(1);
	}
	

	vector<complex<double>> d;					//direction of propagation of the plane wave Init: {0.5, 0}
	if (args["s"].nargs() == 2) {
		d.push_back({ (double)args["s"].as_float(0), 0 });
		d.push_back({ (double)args["s"].as_float(1), 0 });
	}

	double l0 = (double)args["l"].as_float(0);		//wavelength.Init: l0 = 5;
	double k0 = 2 * PI / l0;						//Calculate the free-space wavenumber.

	complex<double> Ex = { (double)args["Ex"].as_float(0), (double)args["Ex"].as_float(1) };	//Input E. Init: {1, 1, 0}
	complex<double> Ey = { (double)args["Ey"].as_float(0), (double)args["Ey"].as_float(1) };
	vector<complex<double>> E0;
	E0.push_back(Ex);
	E0.push_back(Ey);

	//const int LAYERS = args["L"].as_int(); //LAYERS Init: 4

	vector<double> depths;
	vector<complex<double>> ns;
	size_t i = 0;
	while(args["n"].is_set() && args["n"].as_float(i)!=0) {											//n is the real part.
		if (args["kappa"].is_set())																	//kappa is the imaginary part.
			ns.push_back({ (double)args["n"].as_float(i), (double)args["kappa"].as_float(i) });		
		else
			ns.push_back({ (double)args["n"].as_float(i), 0 });										//ns <- n + i * kappa
		i++;
	}
	for (size_t j = 0; j < i; j++)
		if (args["z"].is_set())
			depths.push_back((double)args["z"].as_float(i));
		else {
			depths.push_back((double)-100 + j * 200 / i);
		}
	/*---------------------------------------example_layer--------------------------------------------*/

	calculate_layer(outName, &ns, &depths, &E0, &d, k0);

	return 0;
	std::cin.get();
}