/*
 * OMPOP.c
 *
 *  Created on: 2016年4月5日
 *      Author: asus
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixop.h"

double *OMP(double* colum, double** r, int height, int width, int N) {
	double* hat_y=NULL;
	double*r_n = (double*) malloc(height * sizeof(double));
	double* product = createvector(N);
	double** Aug_t = NULL, **aug_y = NULL;
	double** temp = NULL, **temp1 = NULL, **temp3 = NULL;
	double** pos_array = NULL;
	double** column=(double**)malloc(height*sizeof(double*));
	int M = height, i;
	for(i=0;i<height;i++)
		column[i]=(double*)malloc(sizeof(double));
	for(i=0;i<height;i++)
		column[i][0]=colum[i];
	int times, col, pos;
//	int count=0;//Aug_t的长度为times
	for (i = 0; i < height; i++)
		r_n[i] = column[i][0];
	for (times = 0; times < M; times++) {
		for (col = 0; col < N; col++) {
			product[col] = fabs(
					innerproduct(getrow(r, col, height), column, height));
		}
		pos = max(product, N);
		Aug_t = expandmatrix(Aug_t, getcolumn(r, height, pos), height, times);
		formatcolumn(r, height, pos, 0);
		//aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;
		temp = matrixmultiplic(Aug_t, times, matrixturn(Aug_t, height, times),
				times);
		temp1 = inversematrix(temp, times, getA(temp, times));
		temp = matrixturn(Aug_t, height, times);
		temp3 = matrixmultiplic(temp, times, temp1, times);
		aug_y = matrixmultiplic(column, 1, temp3, times);	//aug_y好多行一列
		//r_n=s-Aug_t*aug_y;
		temp = matrixmultiplic(aug_y, 1, Aug_t, height);
		r_n = vectorsub(column, temp, height);
		//pos_array(times)=pos;
		temp = vectoradd(pos_array, times, pos);
		pos_array = temp;	//pos_array是一行好多列
		//if (abs(aug_y(end))^2/norm(aug_y)<0.05)
		if(pow(fabs(aug_y[times][0]),2)/norm(aug_y,times)<0.05) break;
	}
	hat_y=finish(pos_array,aug_y,times,N);
	return hat_y;
}
