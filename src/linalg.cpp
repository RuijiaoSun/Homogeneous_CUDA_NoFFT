#include "linalg.h"

// This file contains a set of wrapper functions that are linked to the corresponding functions in CLAPACK.
extern "C" {
#include "f2c.h"
#include "clapack.h"
#include "cblas.h"
}


void LINALG_zgetrf(
	int M,
	int N,
	std::complex<double>* A,
	int LDA,
	int* IPIV)
{
	integer INFO;
	zgetrf_((integer*)&M, (integer*)&N, (doublecomplex*)A, (integer*)&LDA, (integer*)IPIV, &INFO);
}

void LINALG_zgetri(
	size_t N,
	std::complex<double>* A,
	int LDA,
	int* IPIV)
{
	integer LWORK = -1;
	std::complex<double> WORK[1];
	integer INFO;
	zgetri_((integer*)&N, (doublecomplex*)A, (integer*)&LDA, (integer*)IPIV, (doublecomplex*)WORK, &LWORK, &INFO);
}
void LINALG_inverse(std::complex<double>* A, int N)
{
	int* IPIV = new int[N + (size_t)1];
	integer LWORK = N * N;
	std::complex<double>* WORK = new std::complex<double>[LWORK];
	integer INFO;

	zgetrf_((integer*)&N, (integer*)&N, (doublecomplex*)A, (integer*)&N, (integer*)IPIV, &INFO);
	zgetri_((integer*)&N, (doublecomplex*)A, (integer*)&N, (integer*)IPIV, (doublecomplex*)WORK, &LWORK, &INFO);

	delete[] IPIV;
	delete[] WORK;
}

void LINALG_zgemm(
	const int M,	//A(M*K) B(K*N)
	const int N,
	const int K,
	std::complex<double>* A,
	const int LDA,   //=K
	std::complex<double>* B,
	const int LDB, //=N
	std::complex<double>* C,
	const int LDC) //=columns of C.
{
	std::complex<double> alpha = 1;
	std::complex<double> beta = 0;
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (OPENBLAS_CONST blasint)M, (OPENBLAS_CONST blasint)N, (OPENBLAS_CONST blasint)K,
		&alpha, A, (OPENBLAS_CONST blasint)LDA, B, (OPENBLAS_CONST blasint)LDB, &beta, C, (OPENBLAS_CONST blasint)LDC);

}