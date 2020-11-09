// This file contains a set of wrapper functions that are linked to the corresponding functions in CLAPACK
#include <complex>

//Solve matrix inverse.
void LINALG_inverse(std::complex<double>* A, int N);

//Solve matrix multiplication. C = A * B.
void LINALG_zgemm(
	const int M,	//A(M*K) B(K*N)
	const int N,
	const int K,
	std::complex<double>* A,
	const int LDA,   //=K
	std::complex<double>* B,
	const int LDB, //=N
	std::complex<double>* C,
	const int LDC); //=columns of C.