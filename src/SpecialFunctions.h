/*
 * Copyright (C) 2014 Quanli Wang
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 */ 

#pragma once
#include <vector>
using namespace std;

class SpecialFunctions
{
public:
	SpecialFunctions(void);
	~SpecialFunctions(void);

	static double gammaln(double x);
	static double betaln(double x, double y);

	static double norminv(double p);			//inverse normal cdf
	static double normcdf(double u);

	static double gammainc(double x, double a);
	static double gammacdf(double x,double a, double b);
	static double gammainv(double p,double a, double b);
	static double gammapdf(double x,double a, double b);
	
	static void cmpower2(int nSize, double *px, double* py, double* pResult);
	static void cmrand(int nSize, MTRand& mt, double* pResult);
	static bool gammarand(double a, double b, int nSize, MTRand& mt, vector<double>& result);
	static double gammarand(double a, double b, MTRand& mt);
	static double chi2rand(double a, MTRand& mt);
	static bool betarand(double a, double b, int nSize, MTRand& mt, vector<double>& result);
	static double betarand(double a, double b, MTRand& mt);
	double betapdf(double x, double a, double b,int logspace);
	static unsigned int binorand(int n, double p, MTRand& mt);

	static double gammarand_int(unsigned int a,MTRand& mt);
	static unsigned int poissonrand(double mu, MTRand& mt);
	static unsigned int negative_binomial_rand(double p, double n, MTRand& mt);


	static double log_sum(double a, double b);
	static double log_gamma_rand(double shape, MTRand& mt);
	static void multinomialrand (unsigned int K, unsigned int N, double *p, unsigned int *n, MTRand& mt);
	static int discreterand(int K, double *p,MTRand& mt);
	static int discreterand_norm(int K, double *p, double norm, MTRand& mt);
};
