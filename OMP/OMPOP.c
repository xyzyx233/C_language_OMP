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

#define R 190
#define X 256

double A[X][X];
double L[X+1][X+1];
double U[X][X];
double u[X][X];
double r[X][X];

double Aug_t[R][X];
double Aug_tt[X][X];
double T[190][X];
double temp1[X][X];//Aug_t'*Aug_t
double temp2[X][X];//(Aug_t'*Aug_t)^(-1)
double temp3[R][X];//(Aug_t'*Aug_t)^(-1)*Aug_t'

void testoutn(double** x, int height, int width, char*str)
{
    int i, j;
    FILE *fto = fopen(str, "wb");
    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            fprintf(fto, "%f\t", x[i][j]);
        }
        fprintf(fto, "\n");
    }
    fclose(fto);
}
void testoutnn(double* x, int height, char*str)
{
    int i;
    FILE *fto = fopen(str, "wb");
    for (i = 0; i < height; i++)
    {
        fprintf(fto, "%f\t", x[i]);
    }
    fclose(fto);
}

int maxxx(double p[]) {
	int imax = 0, i = 0;
	for (i = 0;i<X;i++)
		if (p[imax]<p[i])
			imax = i;
	return imax;
}
double absmax(int col, double t[][X], double r_n[]) {
	int i;
	double result = 0;
	for (i = 0;i < R;i++) {
		result += t[i][col] * r_n[i];
	}
	if (result < 0)
		return -result;
	else
		return result;
}
void expm(int pos, double Aug_t[][X], double T[][X],int times) {
	int i;
		for (i = 0;i < R;i++)
			Aug_t[i][times-1] = T[i][pos];
}
void turnm(double tt[][X],int times,double t[][X]) {
	int i = 0,j=0;
	for (i = 0;i < times;i++)
		for (j = 0;j < R;j++)
			tt[i][j]=t[j][i];
}
void getaug_y(double temp3[][X],double *aug_y,double s[],int times) {
	int i, k;
	double sum;
	for (i = 0;i < times;i++) {
		sum = 0;
		for (k = 0;k < R;k++)
			sum += temp3[i][k] * s[k];
		aug_y[i] = sum;
	}
}
void setr_n(double* r_n,double* s, double* aug_y,double Aug_t[][X],int times) {
	double temp[R];
	double sum = 0;
	int i, j;
	for (i = 0;i < R;i++) {
		sum = 0;
		for (j = 0;j < times;j++) {
			sum += Aug_t[i][j]*aug_y[j];
		}
		temp[i] = sum;
	}
	for (i = 0;i < R;i++) {
		r_n[i] = s[i] - temp[i];
	}
}
void matrixmultiplicx(double matrix1[][X], int width, double matrix2[][X],
	int height, int wid2, double ans[][X]) {
	int i, j, k;
	//	矩阵乘法
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			ans[i][j] = 0;
			for (k = 0; k < wid2; k++) {
				ans[i][j] += matrix2[i][k] * matrix1[k][j];
			}
		}
	}
}
void inversematrixx(double A_[][X], int n, double out[][X]) { //计算每一行每一列的每个元素所对应的余子式，组成A*

	//double **out;
	int i, j, k, x;
	double s;
	//out = (double**)malloc(n*sizeof(double*));
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			A[i][j] = A_[i][j];
			L[i][j] = 0;
		}
	}
	for (i = 0; i < n; i++)
	{
		//out[i] = (double*)malloc(n*sizeof(double));
		L[i][i] = 1;
	}
	for (i = 0;i<n;i++)
		for (j = 0;j < n;j++) {
			U[i][j] = 0;
			L[i][j] = 0;
			u[i][j] = 0;
			r[i][j] = 0;
		}
	//直接三角分解
	//首先计算矩阵L的第1列，矩阵U的第1行
	for (i = 0; i < n; i++)
	{
		U[0][i] = A[0][i];
		if (i != 0)
		{
			L[i][0] = A[i][0] / U[0][0];
		}
	}
	//递推求矩阵L第r行和矩阵U第r列
	for (x = 1; x < n; x++)
	{
		for (int i = x; i < n; i++)
		{
			double tmp1 = 0;
			double tmp2 = 0;
			for (int k = 0; k < x; k++)
			{
				tmp1 += L[x][k] * U[k][i];
				if (i > x)
					tmp2 += L[i][k] * U[k][x];
			}
			U[x][i] = A[x][i] - tmp1;
			if (i > x)
				L[i][x] = (A[i][x] - tmp2) / U[x][x];
		}
	}
	/////////////////////求L和U矩阵的逆//////////////////////////////////////////
	for (i = 0; i<n; i++) /*求矩阵U的逆 */
	{
		u[i][i] = 1 / U[i][i];//对角元素的值，直接取倒数
		for (k = i - 1; k >= 0; k--)
		{
			s = 0;
			for (j = k + 1; j <= i; j++)
				s = s + U[k][j] * u[j][i];
			u[k][i] = -s / U[k][k];//迭代计算，按列倒序依次得到每一个值，
		}
	}
	for (i = 0; i<n; i++) //求矩阵L的逆
	{
		r[i][i] = 1; //对角元素的值，直接取倒数，这里为1
		for (k = i + 1; k<n; k++)
		{
			for (j = i; j <= k - 1; j++)
				r[k][i] = r[k][i] - L[k][j] * r[j][i];   //迭代计算，按列顺序依次得到每一个值
		}
	}
	//////////将r和u相乘，得到逆矩阵
	for (i = 0; i<n; i++)
	{
		for (j = 0; j<n; j++)
		{
			out[i][j] = 0;
		}
	}
	for (i = 0; i<n; i++)
	{
		for (j = 0; j<n; j++)
		{
			for (k = 0; k<n; k++)
			{
				out[i][j] += u[i][k] * r[k][j];
			}
		}
	}
}

double *OMP(double* ss, double** rr, int height, int width, int N)
{
	double* aug_y = NULL;
	double s[R];
	double product[X];
	double r_n[R];
	aug_y = (double*)malloc(R * sizeof(double));
	int pos_array[R];
	double* hat_y = (double*)malloc(X * sizeof(double));

	int c = -1;
	double judge1, judge2;
	int i, j;

	int h = height;
	int w = width;
	int n = N;
	int times;
	int pos = 0;
	int p = 1;

	int col;

	for (i = 0;i < R;i++)
		product[i] = 0;
	for (i = 0;i < R;i++) {
		s[i] = ss[i];
		r_n[i] = s[i];
	}
	for (i = 0;i < R;i++) {
		for (j = 0;j < X;j++) {
			T[i][j] = rr[i][j];
		}
	}
	for (i = 0;i < X;i++)
		hat_y[i] = 0;

	for (times = 1;times <= R;times++) {
		for (col = 0;col <n;col++) {
			product[col] = absmax(col, T, r_n);
		}
		pos = maxxx(product);
		expm(pos, Aug_t, T, times);

		for (i = 0;i < R;i++)
			T[i][pos] = 0;

		turnm(Aug_tt, times, Aug_t);
		matrixmultiplicx(Aug_t, times, Aug_tt, times, R, temp1);
		inversematrixx(temp1, times, temp2);
		matrixmultiplicx(Aug_tt, R, temp2, times, times, temp3);
		getaug_y(temp3, aug_y, s, times);

		setr_n(r_n, s, aug_y, Aug_t, times);
		pos_array[times - 1] = pos;

		//判断条件
		judge1 = fabs(aug_y[times - 1])*fabs(aug_y[times - 1]);
		judge2 = judge1 / norm(aug_y, times);
		if (judge2<0.05) {
			c = times;
			break;
		}
	}

	for (i = 0;i < c;i++) {
		hat_y[pos_array[i]] = aug_y[i];
	}

	return hat_y;
}
