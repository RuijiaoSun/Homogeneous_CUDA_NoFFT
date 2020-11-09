#include "layer.h"
#include "linalg.h"   //LAPACKE support for Visual Studio

#include <cusparse.h>
#include <cuda_runtime.h>
//#include "cublas_v2.h"
#include "cusolverSp.h"
/*----------------------------GPU-----------------------------*/


//Cross product.c is result.
void crossProduct(vector<complex<double>>* a, vector<complex<double>>* b,		//The given matrices.
					vector<complex<double>>* c) {								//Matrix to be gotten.
	c->push_back((*a)[1] * (*b)[2] - (*a)[2] * (*b)[1]);
	c->push_back((*a)[2] * (*b)[0] - (*a)[0] * (*b)[2]);
	c->push_back((*a)[0] * (*b)[1] - (*a)[1] * (*b)[0]);
}

//Calculate the norm for a matrix.
complex<double> Norm(vector<complex<double>>* E) {
	complex<double> sum = 0;
	for (unsigned int i = 0; i < E->size(); i++) {
		sum += (*E)[i] * (*E)[i];
	}
	return sqrt(sum);
}

//Normalize matrix.
void Normalize(vector<complex<double>>* E) {
	complex<double> sum = 0;
	for (unsigned int i = 0; i < E->size(); i++) {
		sum += (*E)[i] * (*E)[i];
	}
	for (unsigned int i = 0; i < E->size(); i++) {
		(*E)[i] = (*E)[i] / sqrt(sum);

	}
}

//Orthogonalization.
void orthogonalize(vector<complex<double>>* E_0rtho, vector<complex<double>>* E0, vector<complex<double>>* d) {
	vector<complex<double>> s;
	if (d->size() == 2) {
		complex<double> dz = sqrt(1.0+0i - pow((*d)[0], 2) - pow((*d)[1], 2));
		d->push_back(dz);
	}
	crossProduct(E0, d, &s);
	crossProduct(d, &s, E_0rtho);
	vector<complex<double>>().swap(s);
}

/*--------------------------------------------------Define Class layersample.--------------------------------------------------------*/

/*Do not try to replace "int" as "size_t". 
This will result in a bunch of warnings and if we continuously change the type of M_rowInd and M_colInd, the EXCEPTION will occur again.*/
size_t layersample::ii(size_t l, size_t c, size_t d) {		//ii(l, c, d) means the column indexes for every element.
	return l * 6 + d * 3 + c - 3;
}				

void layersample::generate_linsys(size_t LAYERS,
									vector<complex<double>>& M,				//All non-zero values in "A" matirx.(A * X = b)
									vector<complex<double>>& b,				//The right hand side column vector.
									vector<complex<double>>& E,				//orthogonalized E0 vectors
									vector<complex<double>>* P,
									bool CPU_op) {			//Solution of the matrices multiplication.
	//Calculate the sz component for each layer.
	s.clear();										//s is the plane wave direction scaled by the refractive index.
	for (size_t i = 0; i < 2; i++)
		s.push_back(d[i] * n[0]);
	sz.clear();
	for (int l = 0; l < LAYERS; l++) {
		sz.push_back(sqrt(pow(n[l], 2) - pow(s[0], 2) - pow(s[1], 2)));
	}

	if (!CPU_op){
		//Computer in GPU.
		vector<int> M_rowInd;												//Sparse matrix M CSR ->row index
		vector<int> M_colInd;												//Sparse matrix M CSR ->number of elements
		M_rowInd.push_back(0);
		////Build M by setting constraints based on Gauss's Law.
		for (size_t l = 0; l < LAYERS; l++) {
			//Set the upward components for each layer.
			//Layer "LAYERS-1" doesn't have a upward component.
			if (l != LAYERS - 1) {
				M.push_back(s[0]);
				M_colInd.push_back((int)ii(l, 0, 1));
				M.push_back(s[1]);
				M_colInd.push_back((int)ii(l, 1, 1));
				M.push_back(-sz[l]);
				M_colInd.push_back((int)ii(l, 2, 1));
				M_rowInd.push_back((int)M.size());
				b.push_back(0);
			}
			//Set the downward components for each layer.
			if (l != 0) {
				M.push_back(s[0]);
				M_colInd.push_back((int)ii(l, 0, 0));
				M.push_back(s[1]);
				M_colInd.push_back((int)ii(l, 1, 0));
				M.push_back(sz[l]);
				M_colInd.push_back((int)ii(l, 2, 0));
				M_rowInd.push_back((int)M.size());
				b.push_back(0);
			}
		}
		//Continue to build M by enforcing a continuous field across boundaries.
		complex<double> arg, arg_in, B;		
		for (size_t l = 1; l < LAYERS; l++) {
			complex<double> sz0 = sz[l - 1];
			complex<double> sz1 = sz[l];

			//Representation of A = np.exp(1j * k0 * sz0 * (self.z[l] - self.z[l - 1]))
			complex<double> A_in = k * sz0 * (z[l] - z[l - 1]);
			complex<double> A_in2 = { -A_in.imag(), A_in.real() };
			complex<double> A = exp(A_in2);

			if (l < LAYERS - 1) {
				double dl = z[l] - z[l + 1];
				arg_in = -k * sz1 * (complex<double>)dl;
				arg = { -arg_in.imag(), arg_in.real() };
				B = exp(arg);
			}
			//if this is the second layer, use the simplified equations that account for the incident field
			if (l == 1) {
				M.push_back(1);
				M_colInd.push_back((int)ii(0, 0, 1));
				M.push_back(-1);
				M_colInd.push_back((int)ii(1, 0, 0));
				if (LAYERS > 2) {
					M.push_back(-B);
					M_colInd.push_back((int)ii(1, 0, 1));
				}
				M_rowInd.push_back((int)M.size());
				b.push_back(-A * E[0]);

				M.push_back(1);
				M_colInd.push_back((int)ii(0, 1, 1));
				M.push_back(-1);
				M_colInd.push_back((int)ii(1, 1, 0));
				if (LAYERS > 2) {
					M.push_back(-B);
					M_colInd.push_back((int)ii(1, 1, 1));
				}
				M_rowInd.push_back((int)M.size());
				b.push_back(-A * E[l]);

				M.push_back(sz0);
				M_colInd.push_back((int)ii(0, 1, 1));
				M.push_back(s[1]);
				M_colInd.push_back((int)ii(0, 2, 1));
				M.push_back(sz1);
				M_colInd.push_back((int)ii(1, 1, 0));
				M.push_back(-s[1]);
				M_colInd.push_back((int)ii(1, 2, 0));
				if (LAYERS > 2) {
					M.push_back(-B * sz1);
					M_colInd.push_back((int)ii(1, 1, 1));
					M.push_back(-B * s[1]);
					M_colInd.push_back((int)ii(1, 2, 1));
				}
				M_rowInd.push_back((int)M.size());
				b.push_back(A * sz0 * E[1] - A * s[1] * E[2]);

				M.push_back(-sz0);
				M_colInd.push_back((int)ii(0, 0, 1));
				M.push_back(-s[0]);
				M_colInd.push_back((int)ii(0, 2, 1));
				M.push_back(-sz1);
				M_colInd.push_back((int)ii(1, 0, 0));
				M.push_back(s[0]);
				M_colInd.push_back((int)ii(1, 2, 0));
				if (LAYERS > 2) {
					M.push_back(B * sz1);
					M_colInd.push_back((int)ii(1, 0, 1));
					M.push_back(B* s[0]);
					M_colInd.push_back((int)ii(1, 2, 1));
				}
				M_rowInd.push_back((int)M.size());
				b.push_back(A * s[0] * E[2] - A * sz0 * E[0]);
			}
			else if (l == LAYERS - 1) {
				M.push_back(A);
				M_colInd.push_back((int)ii(l - 1, 0, 0));
				M.push_back(1);
				M_colInd.push_back((int)ii(l - 1, 0, 1));
				M.push_back(-1);
				M_colInd.push_back((int)ii(l, 0, 0));
				M_rowInd.push_back((int)M.size());
				b.push_back(0);

				M.push_back(A);
				M_colInd.push_back((int)ii(l - 1, 1, 0));
				M.push_back(1);
				M_colInd.push_back((int)ii(l - 1, 1, 1));
				M.push_back(-1);
				M_colInd.push_back((int)ii(l, 1, 0));
				M_rowInd.push_back((int)M.size());
				b.push_back(0);

				M.push_back(-A * sz0);
				M_colInd.push_back((int)ii(l - 1, 1, 0));
				M.push_back(A * s[1]);
				M_colInd.push_back((int)ii(l - 1, 2, 0));
				M.push_back(sz0);
				M_colInd.push_back((int)ii(l - 1, 1, 1));
				M.push_back(s[1]);
				M_colInd.push_back((int)ii(l - 1, 2, 1));
				M.push_back(sz1);
				M_colInd.push_back((int)ii(l, 1, 0));
				M.push_back(-s[1]);
				M_colInd.push_back((int)ii(l, 2, 0));
				M_rowInd.push_back((int)M.size());
				b.push_back(0);

				M.push_back(A * sz0);
				M_colInd.push_back((int)ii(l - 1, 0, 0));
				M.push_back(-A * s[0]);
				M_colInd.push_back((int)ii(l - 1, 2, 0));
				M.push_back(-sz0);
				M_colInd.push_back((int)ii(l - 1, 0, 1));
				M.push_back(-s[0]);
				M_colInd.push_back((int)ii(l - 1, 2, 1));
				M.push_back(-sz1);
				M_colInd.push_back((int)ii(l, 0, 0));
				M.push_back(s[0]);
				M_colInd.push_back((int)ii(l, 2, 0));
				M_rowInd.push_back((int)M.size());
				b.push_back(0);
			}
			else {
				M.push_back(A);
				M_colInd.push_back((int)ii(l - 1, 0, 0));
				M.push_back(1);
				M_colInd.push_back((int)ii(l - 1, 0, 1));
				M.push_back(-1);
				M_colInd.push_back((int)ii(l, 0, 0));
				M.push_back(-B);
				M_colInd.push_back((int)ii(l, 0, 1));
				M_rowInd.push_back((int)M.size());
				b.push_back(0);

				M.push_back(A);
				M_colInd.push_back((int)ii(l - 1, 1, 0));
				M.push_back(1);
				M_colInd.push_back((int)ii(l - 1, 1, 1));
				M.push_back(-1);
				M_colInd.push_back((int)ii(l, 1, 0));
				M.push_back(-B);
				M_colInd.push_back((int)ii(l, 1, 1));
				M_rowInd.push_back((int)M.size());
				b.push_back(0);

				M.push_back(-A * sz0);
				M_colInd.push_back((int)ii(l - 1, 1, 0));
				M.push_back(A * s[1]);
				M_colInd.push_back((int)ii(l - 1, 2, 0));
				M.push_back(sz0);
				M_colInd.push_back((int)ii(l - 1, 1, 1));
				M.push_back(s[1]);
				M_colInd.push_back((int)ii(l - 1, 2, 1));
				M.push_back(sz1);
				M_colInd.push_back((int)ii(l, 1, 0));
				M.push_back(-s[1]);
				M_colInd.push_back((int)ii(l, 2, 0));
				M.push_back(-B * sz1);
				M_colInd.push_back((int)ii(l, 1, 1));
				M.push_back(-B * s[1]);
				M_colInd.push_back((int)ii(l, 2, 1));
				M_rowInd.push_back((int)M.size());
				b.push_back(0);

				M.push_back(A * sz0);
				M_colInd.push_back((int)ii(l - 1, 0, 0));
				M.push_back(-A * s[0]);
				M_colInd.push_back((int)ii(l - 1, 2, 0));
				M.push_back(-sz0);
				M_colInd.push_back((int)ii(l - 1, 0, 1));
				M.push_back(-s[0]);
				M_colInd.push_back((int)ii(l - 1, 2, 1));
				M.push_back(-sz1);
				M_colInd.push_back((int)ii(l, 0, 0));
				M.push_back(s[0]);
				M_colInd.push_back((int)ii(l, 2, 0));
				M.push_back(B * sz1);
				M_colInd.push_back((int)ii(l, 0, 1));
				M.push_back(B * s[0]);
				M_colInd.push_back((int)ii(l, 2, 1));
				M_rowInd.push_back((int)M.size());
				b.push_back(0);
			}
		}
		cudaError_t cudaStatus;
		cusolverStatus_t cusolverStatus;
		cusparseStatus_t cusparseStatus;
		cusolverSpHandle_t handle = NULL;
		cusparseHandle_t cusparseHandle = NULL;
		cudaStream_t stream = NULL;
		cusparseMatDescr_t descrM = NULL;
		cuDoubleComplex * csrValM_, * b_, *P_;
		size_t rowsA = b.size(), colsA = b.size(), nnA = M.size(), baseM_ = 0;		//nnA is the number of non-zero elements.
		int* csrRowPtrM = NULL;														//row index M_rowInd projected to GPU.
		int* csrColIndM = NULL; //CSR(A) from I/O.									// M_colInd projected to GPU.
		double tol = 1.e-12; int reorder = 0;
		int singularity = 0;
	
		//Initialize.
		cusolverStatus = cusolverSpCreate(&handle);
		int num = 1;
		cudaStatus = cudaGetDevice(&num);
		cusparseStatus = cusparseCreate(&cusparseHandle);
		cudaStatus = cudaStreamCreate(&stream);
		cusolverStatus = cusolverSpSetStream(handle, stream);
		cusparseStatus = cusparseSetStream(cusparseHandle, stream);
		cusparseStatus = cusparseCreateMatDescr(&descrM);
		cusparseStatus = cusparseSetMatType(descrM, CUSPARSE_MATRIX_TYPE_GENERAL);
		if (baseM_) {
			cusparseStatus = cusparseSetMatIndexBase(descrM, CUSPARSE_INDEX_BASE_ONE);
		}
		else {
			cusparseStatus = cusparseSetMatIndexBase(descrM, CUSPARSE_INDEX_BASE_ZERO);
		}

		cudaStatus = cudaMalloc((void**)&csrRowPtrM, sizeof(int) * (rowsA + 1));				//Projection of M_rowInd.
		cudaStatus = cudaMalloc((void**)&csrColIndM, sizeof(int) * M_colInd.size());			//Projection of M_colInd.
		cudaStatus = cudaMalloc((void**)&csrValM_, sizeof(cuDoubleComplex) * M.size());			//Projection of M.
		cudaStatus = cudaMalloc((void**)&b_, sizeof(cuDoubleComplex) * b.size());				//Projection of b.
		cudaStatus = cudaMalloc((void**)&P_, sizeof(cuDoubleComplex) * b.size());				//Projection of P.

		cudaStatus = cudaMemcpy(csrValM_, M.data(), M.size() * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
		cudaStatus = cudaMemcpy(csrRowPtrM, M_rowInd.data(), M_rowInd.size() * sizeof(int), cudaMemcpyHostToDevice);
		cudaStatus = cudaMemcpy(csrColIndM, M_colInd.data(), M_colInd.size() * sizeof(int), cudaMemcpyHostToDevice);
		cudaStatus = cudaMemcpy(b_, b.data(), b.size() * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
		// Output the current CUDA error.
		//if (cudaStatus != cudaSuccess) {
		//	cout<<"%s " << cudaGetErrorString(cudaStatus) << endl;
		//}
		P->resize(rowsA);							//P is the to-be-solved matrix in CPU.
		//QR method.
		cusolverStatus = cusolverSpZcsrlsvqr(handle, (int)rowsA, (int)nnA, descrM, csrValM_, csrRowPtrM, csrColIndM, b_, tol, reorder, P_, (int*)&singularity);
		/*cusparseStatus = cusparseZsctr(*cusparseHandle, rowsA, g_z, g_Q, g_x, CUSPARSE_INDEX_BASE_ZERO);*/
		cudaStatus = cudaMemcpyAsync(P->data(), P_, sizeof(cuDoubleComplex) * rowsA, cudaMemcpyDeviceToHost, stream);
	
		cudaStatus = cudaFree(csrRowPtrM);
		cudaStatus = cudaFree(csrColIndM);
		cudaStatus = cudaFree(csrValM_);
		cudaStatus = cudaFree(b_);
		cudaStatus = cudaFree(P_);
		vector<int>().swap(M_rowInd);
		vector<int>().swap(M_colInd);
	}
	else {
		//Work on CPU.
		M.resize(6 * (LAYERS - 1) * 6 * (LAYERS - 1));
		b.resize(6 * (LAYERS - 1));

		size_t ei = 0;
		//Set constraints based on Gauss's Law.
		for (size_t l = 0; l < LAYERS; l++) {
			//Set the upward components for each layer.
			//Layer "LAYERS-1" doesn't have a upward component.
			if (l != LAYERS - 1) {
				M[ei * 6 * (LAYERS - 1) + ii(l, 0, 1)] = s[0];
				M[ei * 6 * (LAYERS - 1) + ii(l, 1, 1)] = s[1];
				M[ei * 6 * (LAYERS - 1) + ii(l, 2, 1)] = -sz[l];
				ei += 1;
			}
			//Set the downward components for each layer.
			if (l != 0) {
				M[ei * 6 * (LAYERS - 1) + ii(l, 0, 0)] = s[0];
				M[ei * 6 * (LAYERS - 1) + ii(l, 1, 0)] = s[1];
				M[ei * 6 * (LAYERS - 1) + ii(l, 2, 0)] = sz[l];
				ei += 1;
			}
		}
		//Enforce a continuous field across boundaries.
		complex<double> arg, arg_in, B;
		for (size_t l = 1; l < LAYERS; l++) {
			complex<double> sz0 = sz[l - 1];
			complex<double> sz1 = sz[l];

			//Representation of A = np.exp(1j * k0 * sz0 * (self.z[l] - self.z[l - 1]))
			complex<double> A_in = k * sz0 * (z[l] - z[l - 1]);
			complex<double> A_in2 = { -A_in.imag(), A_in.real() };
			complex<double> A = exp(A_in2);

			if (l < LAYERS - 1) {
				double dl = z[l] - z[l + 1];
				arg_in = -k * sz1 * (complex<double>)dl;
				arg = { -arg_in.imag(), arg_in.real() };
				B = exp(arg);
			}
			//if this is the second layer, use the simplified equations that account for the incident field
			if (l == 1) {
				M[ei * 6 * (LAYERS - 1) + ii(0, 0, 1)] = 1;
				M[ei * 6 * (LAYERS - 1) + ii(1, 0, 0)] = -1;
				if (LAYERS > 2) {
					M[ei * 6 * (LAYERS - 1) + ii(1, 0, 1)] = -B;
				}
				b[ei] = -A * E[0];
				ei += 1;

				M[ei * 6 * (LAYERS - 1) + ii(0, 1, 1)] = 1;
				M[ei * 6 * (LAYERS - 1) + ii(1, 1, 0)] = -1;
				if (LAYERS > 2) {
					M[ei * 6 * (LAYERS - 1) + ii(1, 1, 1)] = -B;
				}
				b[ei] = -A * E[l];
				ei += 1;

				M[ei * 6 * (LAYERS - 1) + ii(0, 2, 1)] = s[1];
				M[ei * 6 * (LAYERS - 1) + ii(0, 1, 1)] = sz0;
				M[ei * 6 * (LAYERS - 1) + ii(1, 2, 0)] = -s[1];
				M[ei * 6 * (LAYERS - 1) + ii(1, 1, 0)] = sz1;
				if (LAYERS > 2) {
					M[ei * 6 * (LAYERS - 1) + ii(1, 2, 1)] = -B * s[1];
					M[ei * 6 * (LAYERS - 1) + ii(1, 1, 1)] = -B * sz1;
				}
				b[ei] = A * sz0 * E[1] - A * s[1] * E[2];
				ei += 1;

				M[ei * 6 * (LAYERS - 1) + ii(0, 0, 1)] = -sz0;
				M[ei * 6 * (LAYERS - 1) + ii(0, 2, 1)] = -s[0];
				M[ei * 6 * (LAYERS - 1) + ii(1, 0, 0)] = -sz1;
				M[ei * 6 * (LAYERS - 1) + ii(1, 2, 0)] = s[0];
				if (LAYERS > 2) {
					M[ei * 6 * (LAYERS - 1) + ii(1, 0, 1)] = B * sz1;
					M[ei * 6 * (LAYERS - 1) + ii(1, 2, 1)] = B * s[0];
				}
				b[ei] = A * s[0] * E[2] - A * sz0 * E[0];
				ei += 1;
			}
			else if (l == LAYERS - 1) {
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 0, 0)] = A;
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 0, 1)] = 1;
				M[ei * 6 * (LAYERS - 1) + ii(l, 0, 0)] = -1;
				ei += 1;

				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 1, 0)] = A;
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 1, 1)] = 1;
				M[ei * 6 * (LAYERS - 1) + ii(l, 1, 0)] = -1;
				ei += 1;

				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 2, 0)] = A * s[1];
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 1, 0)] = -A * sz0;
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 2, 1)] = s[1];
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 1, 1)] = sz0;
				M[ei * 6 * (LAYERS - 1) + ii(l, 2, 0)] = -s[1];
				M[ei * 6 * (LAYERS - 1) + ii(l, 1, 0)] = sz1;
				ei += 1;

				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 0, 0)] = A * sz0;
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 2, 0)] = -A * s[0];
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 0, 1)] = -sz0;
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 2, 1)] = -s[0];
				M[ei * 6 * (LAYERS - 1) + ii(l, 0, 0)] = -sz1;
				M[ei * 6 * (LAYERS - 1) + ii(l, 2, 0)] = s[0];
				ei += 1;
			}
			else {
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 0, 0)] = A;
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 0, 1)] = 1;
				M[ei * 6 * (LAYERS - 1) + ii(l, 0, 0)] = -1;
				M[ei * 6 * (LAYERS - 1) + ii(l, 0, 1)] = -B;
				ei += 1;

				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 1, 0)] = A;
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 1, 1)] = 1;
				M[ei * 6 * (LAYERS - 1) + ii(l, 1, 0)] = -1;
				M[ei * 6 * (LAYERS - 1) + ii(l, 1, 1)] = -B;
				ei += 1;

				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 2, 0)] = A * s[1];
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 1, 0)] = -A * sz0;
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 2, 1)] = s[1];
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 1, 1)] = sz0;
				M[ei * 6 * (LAYERS - 1) + ii(l, 2, 0)] = -s[1];
				M[ei * 6 * (LAYERS - 1) + ii(l, 1, 0)] = sz1;
				M[ei * 6 * (LAYERS - 1) + ii(l, 2, 1)] = -B * s[1];
				M[ei * 6 * (LAYERS - 1) + ii(l, 1, 1)] = -B * sz1;
				ei += 1;

				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 0, 0)] = A * sz0;
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 2, 0)] = -A * s[0];
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 0, 1)] = -sz0;
				M[ei * 6 * (LAYERS - 1) + ii(l - 1, 2, 1)] = -s[0];
				M[ei * 6 * (LAYERS - 1) + ii(l, 0, 0)] = -sz1;
				M[ei * 6 * (LAYERS - 1) + ii(l, 2, 0)] = s[0];
				M[ei * 6 * (LAYERS - 1) + ii(l, 0, 1)] = B * sz1;
				M[ei * 6 * (LAYERS - 1) + ii(l, 2, 1)] = B * s[0];
				ei += 1;
			}
		}

		complex<double>* M_ = new complex<double>[M.size()];
		complex<double>* b_ = new complex<double>[b.size()];
		complex<double>* P_ = new complex<double>[b.size()];
		for (size_t i = 0; i < M.size(); i++) {
			M_[i] = M[i];
			if (i < b.size())	b_[i] = b[i];
		}
		LINALG_inverse(M_, (int)(6 * (LAYERS - 1)));
		LINALG_zgemm((int)(6 * (LAYERS - 1)), (int)1, (int)(6 * (LAYERS - 1)), M_, (int)(6 * (LAYERS - 1)), b_, (int)1, P_, (int)1);
		for (int i = 0; i < b.size(); i++) {
			P->push_back(P_[i]);
		}

		delete[] M_;
		delete[] b_;
		delete[] P_;
	}
}


//Build matrix and get E.
void layersample::solve(vector<complex<double>>* E, bool CPU_op) {		//orthogonalized E0 vectors.
	size_t LAYERS = n.size();
	//Store the matrix and RHS vector.
	vector<complex<double>> M;						//All non-zero values in the sparse matrix.
	vector<complex<double>> b;						//The right hand side column vector.

	//Evaluate the linear system.
	vector<complex<double>> P;						//Solution of matrix.
	layersample::generate_linsys(LAYERS, M, b, *E, &P, CPU_op);

	//Store the coefficients for each layer.
	//Pt[3, L] transmission. Pr[3, L] reflection.
	Pt.resize(3 * LAYERS); 
	Pr.resize(3 * LAYERS);

	for (size_t l = 0; l < LAYERS; l++) {
		if (l == 0) {
			Pt[0] = (complex<double>)(*E)[0];
			Pt[LAYERS] = (complex<double>)(*E)[1];
			Pt[2 * LAYERS] = (complex<double>)(*E)[2];
		}
		else {
			Pt[l] = P[ii(l, 0, 0)];
			Pt[l + LAYERS] = P[ii(l, 1, 0)];
			Pt[l + 2 * LAYERS] = P[ii(l, 2, 0)];
		}

		if (l == LAYERS - 1) {
			Pr[LAYERS - 1] = 0;
			Pr[2 * LAYERS - 1] = 0;
			Pr[3 * LAYERS - 1] = 0;
		}
		else {
			Pr[l] = P[ii(l, 0, 1)];
			Pr[l + LAYERS] = P[ii(l, 1, 1)];
			Pr[l + 2 * LAYERS] = P[ii(l, 2, 1)];
		}
	}
	vector<complex<double>>().swap(M);
	vector<complex<double>>().swap(b);
	vector<complex<double>>().swap(P);
}