/*
 * matrixop.c
 *
 *  Created on: 2016年4月5日
 *      Author: asus
 */

#include <math.h>
#include  <stdio.h>
#include <stdlib.h>

double** matrixmultiplic(double** matrix1, int width, double** matrix2,
		int height) {
	int i, j, k;
	double** ans = NULL;
	//double** ans = creatmatrix(height,width); 这种写法会报错 conflicting types for 原因未知(目前)
	ans = (double**) malloc(height * sizeof(double*));
	if (ans == NULL)
		return NULL;
	for (i = 0; i < height; i++) {
		ans[i] = (double*) malloc(width * sizeof(double));
		if (ans[i] == NULL)
			return NULL;
	}
//	矩阵乘法
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			ans[i][j] = 0;
			for (k = 0; k < width; k++) {
				ans[i][j] += matrix2[i][k] * matrix1[k][j];
			}
		}
	}
	return ans;
}
void clear(double** matrix, int matirxsize) {
	int i;
	for (i = 0; i < matirxsize; i++) {
		free(matrix[i]);
	}
	free(matrix);
}
double** matrixturn(double** matrix, int height, int width) {
	int i = 0, j;
	double** turn = (double**) malloc(width * sizeof(double*));
	if (!turn)
		return NULL;
	for (i = 0; i < width; i++)
		turn[i] = (double*) malloc(height * sizeof(double));
	for (i = 0; i < width; i++)
		for (j = 0; j < height; j++)
			turn[i][j] = matrix[j][i];
	return turn;
}
double** creatmatrix(int height, int width) {
	double **pixel = (double**) malloc(height * sizeof(double*));
	if (pixel == NULL) {
		printf("It is out of memory!\n");
		return NULL;
	}
	int i, j;
	for (i = 0; i < height; i++) {
		pixel[i] = (double*) malloc(width * sizeof(double));
		if (pixel[i] == NULL) {
			printf("It is out of memory!\n");
			return NULL;
		}
	}
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			pixel[i][j] = 0;
		}
	}
	return pixel;
}
double* getcolumn(double** matrix, int height, int col) {
	double* ans = NULL;
	ans = (double*) malloc(height * sizeof(double));
	if (!ans)
		return NULL;
	int i;
	for (i = 0; i < height; i++) {
		ans[i] = matrix[i][col];
	}
	return ans;
}
void setcolum(double** matrix, int column, double* vactor, int row) {
	int i;
	for (i = 0; i < row; i++) {
		matrix[column][i] = vactor[i];
	}
}
double* createvector(int n) {
	double* ans = (double*) malloc(n * sizeof(double));
	int i = 0;
	if (!ans)
		return NULL;
	for (i = 0; i < n; i++)
		ans[i] = 0;
	return ans;
}
double innerproduct(double* a, double** b, int len) {
	int i;
	double sum = 0;
	for (i = 0; i < len; i++) {
		sum += a[i] * b[i][0];
	}
	return sum;
}
double* getrow(double** matrix, int column, int row) {
	int i;
	double* vactor = createvector(column);
	for (i = 0; i < row; i++) {
		matrix[i][column] = vactor[i];
	}
	return vactor;
}
int max(double* vector, int num) {
	int i, max = 0;
	for (i = 1; i < num; i++) {
		if (vector[i] > vector[i - 1]) {
			max = i;
		}
	}
	return max;
}
void formatcolumn(double** matrix, int height, int column, int val) {
	int i;
	for (i = 0; i < height; i++) {
		matrix[i][column] = val;
	}
}
double** expandmatrix(double **ma1, double *ve, int height, int width) {
	int i, j;
	double** ma2 = (double**) malloc(height * sizeof(double*));
	for (i = 0; i < height; i++)
		ma2[i] = (double*) malloc((width + 1) * sizeof(double));
	for (i = 0; i < height; i++)
		for (j = 0; j < width; j++)
			ma2[i][j] = ma1[i][j];
	for (i = 0; i < height; i++) {
		ma2[i][width] = ve[i];
	}
	return ma2;
}
double getA(double** arcs, int n) //按第一行展开计算|A|
{
	if (n == 1) {
		return arcs[0][0];
	}
	int i, j, k;
	double det = 0;
	double** temp = (double**) malloc(n * sizeof(double*));
	for (i = 0; i < n; i++) {
		temp[i] = (double*) malloc(n * sizeof(double));
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n - 1; j++) {
			for (k = 0; k < n - 1; k++) {
				temp[j][k] = arcs[j + 1][(k >= i) ? k + 1 : k];
			}
		}
		double t = getA(temp, n - 1);
		if (i % 2 == 0) {
			det += arcs[0][i] * t;
		} else {
			det -= arcs[0][i] * t;
		}
	}
	return det;
}
double** inversematrix(double** arcs, int n, double det) { //计算每一行每一列的每个元素所对应的余子式，组成A*
	double** ans = NULL;
	int i, j, k, t;
	if (n == 1) {
		ans = (double**) malloc(sizeof(double*));
		ans[0] = (double*) malloc(sizeof(double));
		ans[0][0] = 1;
		return ans;
	}
	double **temp = (double**) malloc(n * sizeof(double*));
	ans = (double**) malloc(sizeof(n * sizeof(double*)));
	for (i = 0; i < n; i++) {
		temp[i] = (double*) malloc(n * sizeof(double));
		ans[i] = (double*) malloc(n * sizeof(double));
	}
	//求矩阵的逆
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			ans[i][j] = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n - 1; k++) {
				for (t = 0; t < n - 1; t++) {
					temp[k][t] = arcs[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
				}
			}
			ans[j][i] = getA(temp, n - 1);
			if ((i + j) & 1) {
				ans[j][i] = -ans[j][i];
			}
		}
	}
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			ans[i][j] = ans[i][j] / det;
	return ans;
}
double* vectorsub(double** vec1, double** vec2, int num) {
	double* ans = (double*) malloc(num * sizeof(double));
	int i;
	for (i = 0; i < num; i++)
		ans[i] = vec1[i][0] - vec2[0][i];
	return ans;
}
double** vectoradd(double** vec, int num, double n) {
	double** ans=(double**)malloc(sizeof(double*));
	ans[0]= (double*) realloc(ans[0], (num + 1) * sizeof(double));
	ans[0][num] = n;
	return ans;
}
double norm(double* vec,int num){
	double ans;
	int i;
	for(i=0;i<num;i++)
		ans+=vec[i]*vec[i];
	return sqrt(ans);
}
double* finish(double** pos, double** aug, int num,int n){
	double* ans=createvector(n);
	int i;
	for(i=0;i<num;i++){
		ans[(int)pos[0][i]]=aug[i][0];
	}
	return ans;
}
