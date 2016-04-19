/*
 * matrixop.h
 *
 *  Created on: 2016年4月5日
 *      Author: eric
 */

#ifndef MATRIXOP_H_
#define MATRIXOP_H_

double** matrixmultiplic(double** matrix1, int width, double** matrix2,
		int height,int wid2);	//后边的是矩阵相乘中前边的矩阵
void clear(double** matrix, int matirxsize);
double** matrixturn(double** matrix, int height, int width);
double** creatmatrix(int height, int width);	//生成指定大小的矩阵
double* createvector(int len);
double* getcolumn(double** matrix,int height, int col);
void setcolum(double** matrix, int column,double* vactor,int row);
void formatcolumn(double** matrix, int height, int column, double val);
double* getrow(double** matrix, int row,int column);
double innerproduct(double* a,double** b, int len);
int max(double* vector, int num);
double** expandmatrix(double **ma1, double *ve, int height, int width);
double** inversematrix(double** arcs, int n, double det);
double getA(double** arcs, int n);
double* vectorsub(double** vec1, double** vec2, int num);
double** vectoradd(double** vec,int num,double n);
double norm(double** vec,int num);
double* finish(double** pos, double** aug, int num, int n);
void RightShift(double *arr, int N, int k);

#endif /* MATRIXOP_H_ */
