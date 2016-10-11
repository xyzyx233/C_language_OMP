/*
 * matrixop.h
 *
 *  Created on: 2016��4��5��
 *      Author: eric
 */

#ifndef MATRIXOP_H_
#define MATRIXOP_H_

double** matrixmultiplic(double** matrix1, int width, double** matrix2,
		int height,int wid2);	//��ߵ��Ǿ��������ǰ�ߵľ���
void clear(double** matrix, int matirxsize);
double** matrixturn(double** matrix, int height, int width);
double** creatmatrix(int height, int width);	//����ָ����С�ľ���
double* createvector(int len);
double* getcolumn(double** matrix,int height, int col);
void setcolum(double** matrix, int column,double* vactor,int row);
void formatcolumn(double** matrix, int height, int column, double val);
double* getrow(double** matrix, int row,int column);
double innerproduct(double* a,double* b, int len);
int maxx(double* vector, int num);
void expandmatrix(double **ma1, double *ve, int height, int width);
double** inversematrix(double **A_, int n);
double getA(double** arcs, int n);
double* vectorsub(double** vec1, double** vec2, int num);
void vectoradd(int* vec,int num,int n);
double norm(double* vec,int num);
double* finish(int* pos, double** aug, int num, int n);
void RightShift(double *arr, int N, int k);

#endif /* MATRIXOP_H_ */