/*
Copyright (c) 2010, Matt Berger
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list
  of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list
  of conditions and the following disclaimer in the documentation and/or other
  materials provided with the distribution.
- Neither the name of the University of Utah nor the names of its contributors may be used
  to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef DENSEEIGENSYSTEM_H
#define DENSEEIGENSYSTEM_H

#include "DenseMatrix.h"

#define CHECK_EIGENSYSTEM false

class DenseEigenSystem  {
	public:
		DenseEigenSystem(DenseMatrix *_matrix, int _il, int _iu);
		~DenseEigenSystem();

		BLinAlg::VectorBLA** getEigenvectors()  {
			return eigenvectors;
		}

		BLinAlg::VectorBLA* getEigenvector(int _i)  {
			return eigenvectors[_i-il+1];
		}

		double* getEigenvalues()  {
			return eigenvalues;
		}

		double getEigenvalue(int _i)  {
			return eigenvalues[_i-il+1];
		}

		int numEigens()  {
			return systemSize;
		}

		double determinant()  {
			double determ = 1.0;
			for(int i = 0; i < systemSize; i++)
				determ *= eigenvalues[i];
			return determ;
		}

		void solve();

	private:
		DenseMatrix* matrix;

		int il, iu;
		int systemSize;

		double* eigenvalues;
		BLinAlg::VectorBLA** eigenvectors;
};

#endif
