/*
 * Copyright (C) 2007-2014 Daniel Manrique-Vallier
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 * Modified by Quanli Wang, 2014
 */ 
#ifndef _CPARAM_H
#define _CPARAM_H

#include "CArrayND.h"
#include <string.h>
#include "MersenneTwister.h"
#include "SpecialFunctions.h"

class CParam {
public:
	CParam(int J, int K, int L, int *levelsJ, int *cumLevelsJ, int n, int Nmis_max, int **MCZ, int nZeroMC, double a_alpha, double b_alpha, int** x)
			: a_alpha(a_alpha), b_alpha(b_alpha){
		class_construct(J, K, L, levelsJ,cumLevelsJ, n);
		class_construct(Nmis_max, MCZ, nZeroMC, x);
	}
	CParam(CData* dat, int K, int Nmis_max,  double a_alpha, double b_alpha)
			: a_alpha(a_alpha), b_alpha(b_alpha){ 
				class_construct(dat->J,  K, dat->L, dat->levelsJ, dat->cumLevelsJ,dat->n);	
		class_construct(Nmis_max, dat->ZeroMC_IJ, dat->nZeroMC, dat->x);
	}
	virtual ~CParam(); //Destructor

	int *zI; //(z_i = k) i=1..N_MAX //just n
	double **psiJKL; //(psi_{jkl}) l -> # of levels
	double *nuK;
	int *countK;

	//structure
	int J, K, L, *levelsJ, n, *cumLevelsJ;

	//auxiliary
	int** aux_dirCumJK;

	//DP parameters 
	double *log_nuK;
	double alpha;
	int k_star;
	//priors
	double a_alpha, b_alpha; // alpha ~ Gamma[a_alpha, b_alpha]

	//partition
	int **xIJ; //for holding data and augmented data.
	int **MCZ; //for holing a local copy of marginal conditions. 
	double *pZeroMC_I, prob_zero;
	int *z2_Nmax;
	int **x2_NMax_J;
	unsigned int *count_partition; //imputed number of individuals in each partition
	int Nmis, N_mis_max, nZeroMC;

	virtual void initizalize(MTRand & mt); //Initialize the parameters
	virtual void predict(int* xIJ_flat, double *ret, int I){
		int * cellJ = xIJ_flat;
		double * p_sum = ret;
		for (int i = 0; i < I; ++i, cellJ += J, ++p_sum){
			*p_sum = 0.0;
			for(int k = 0; k < K; ++k){
				double prod = nuK[k];
				for(int j = 0; j < J; ++j){
					int l = cellJ[j];
					if ( l != -1 ) prod *= psiJKL[cumLevelsJ[j]+l][k];
				}
				*p_sum += prod;
			}
		}

		const double renorm = 1.0 / (1.0 - prob_zero);
		for (double* it = ret; it != ret + I; ++it){
			*it *= renorm;
		}
	}

private:
	void class_construct(int J, int K, int L, int *levelsJ, int *cumLevelsJ, int n);
	void class_construct(int Nmis_max, int** MCZ_, int nZeroMC, int **x);	
};

#endif
