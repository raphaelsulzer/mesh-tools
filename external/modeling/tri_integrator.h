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

#ifndef TRI_INTEGRATOR
#define TRI_INTEGRATOR

#include "mesh_types.h"
#include "ImplicitFunction.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

struct SubTri  {
	SubTri()  {}
	SubTri(Vector3RB _origin, Vector3RB _v2, Vector3RB _v3) : o(_origin) , v2(_v2) , v3(_v3) {}
	~SubTri()  {}
	Vector3RB o, v2, v3;
};

struct TriFunction  {
	TriFunction() : f1(0), f2(0), f3(0) {}
	TriFunction(double _f1, double _f2, double _f3) : f1(_f1), f2(_f2), f3(_f3) {}
	double f1, f2, f3;
};

class TriIntegrator  {
	public:
		TriIntegrator(Vector3RB _point, Vector3RB _v1, Vector3RB _v2, Vector3RB _v3, double _eps, double _lambda);
		TriIntegrator(Vector3RB _point, Vector3RB _v1, Vector3RB _v2, Vector3RB _v3, double _eps);
		TriIntegrator(Vector3RB _point, double _eps);
		~TriIntegrator();

		void set_steps(int _numSteps, double* _steps)  {
			num_steps = _numSteps;
			nonuniform_steps = _steps;
		}

		bool useApproximation()  { return use_approximation; }

		double weight_eval(Vector3RB _p1, Vector3RB _p2);
		double weight_eval_sqd(Vector3RB _p1, Vector3RB _p2);
		double weight_eval_cbd(Vector3RB _p1, Vector3RB _p2);

		void tri_integrals();

		double constant_integration();
		double numerical_constant_integration();

		double linear_integration(TriFunction _f);
		double numerical_linear_integration(TriFunction _f);

		double product_integration(TriFunction _f, TriFunction _g);
		double numerical_product_integration(TriFunction _f, TriFunction _g);

		Vector3RB gradient_integration();
		Vector3RB numerical_gradient_integration();

		double inv_quadratic_integral(double _a, double _b, double _c);
		double inv_sqd_quadratic_integral(double _a, double _b, double _c);
		double inv_linear_sqd_quadratic_integral(double _a, double _b, double _c);

		double inv_cubic_integral(double _a, double _b, double _c);
		double inv_linear_cubic_integral(double _a, double _b, double _c);

		double constant_line_integration(Vector3RB _p1, Vector3RB _p2);
		double linear_line_integration(Vector3RB _p1, Vector3RB _p2, double _f1, double _f2);
		double quadratic_line_integration(Vector3RB _p1, Vector3RB _p2, double _f1, double _f2, double _g1, double _g2);

		Vector3RB gradient_line_integration(Vector3RB _p1, Vector3RB _p2);

		Vector3RB closestPoint()  { return closest_point; }

		bool verify;

	private:
		Vector3RB v1, v2, v3;
		Vector3RB point;

		Vector3RB normal;
		double area;

		double eps;

		int num_subtris;
		SubTri subdivided_tris[3];

		Vector3RB closest_point;
		bool use_approximation;

		int num_steps;
		double* nonuniform_steps;

		Vector3RB ray_plane_intersection(Vector3RB _p, Vector3RB _n, Vector3RB _m, Vector3RB _d);
		Vector3RB orthogonal_projection(Vector3RB _p, Vector3RB _x, Vector3RB _normal);
		double signed_dist(Vector3RB _p, Vector3RB _x, Vector3RB _normal);
		Vector3RB edge_projection(Vector3RB _v1, Vector3RB _v2, Vector3RB _normal, Vector3RB _point);

		void barycentric_coordinate(Vector3RB _point, double& b1, double& b2, double& b3);
};

#endif
