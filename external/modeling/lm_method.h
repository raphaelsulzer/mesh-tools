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

#ifndef LEVENBERGMARQUARDT_H
#define LEVENBERGMARQUARDT_H

#include "DenseMatrix.h"

class LevenbergMarquardt  {
	public:
		LevenbergMarquardt(int _solverIterations, float _energyEps);
		LevenbergMarquardt();
		~LevenbergMarquardt();

		virtual BLinAlg::VectorBLA* initial_guess() = 0;
		virtual BLinAlg::VectorBLA* next_function(BLinAlg::VectorBLA* _approximation) = 0;
		virtual DenseMatrix* next_jacobian(BLinAlg::VectorBLA* _approximation) = 0;

		BLinAlg::VectorBLA* nonlinear_least_squares();

		double squared_energy()  { return sqd_energy; }

	protected:
		int solver_iterations;
		float energy_eps;
		double sqd_energy;
};

#endif
