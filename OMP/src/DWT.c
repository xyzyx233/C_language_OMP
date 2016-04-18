/*
 * DWT.c
 *
 *  Created on: 2016年4月5日
 *      Author: asus
 */

#include "DWT.h"
#include "matrixop.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

double **DWT(int n) {
	int i, j; //外循环控制变量
	int x = 0, y = 0; //内循环控制变量
	double temp1 = 0, temp2 = 0, temp3 = 0, temp4 = 0;	//交换用临时变量
	double h[16] = { -0.0034, -0.0005, 0.0317, 0.0076, -0.1433, -0.0613, 0.4814,
			0.7772, 0.3644, -0.0519, -0.0272, 0.0491, 0.0038, -0.0150, -0.0003,
			0.0019 };	//常量
	double g[16] = { -0.0019, -0.0003, 0.0150, 0.0038, -0.0491, -0.0272, 0.0519,
			0.3644, -0.7772, 0.4814, 0.0613, -0.1433, -0.0076, 0.0317, 0.0005,
			-0.0034 };	//常量
	int l = 16;	//length(h)
	//为结果分配足够的内存空间，如果分配失败则返回NULL，并赋值
	double **ans = (double**) malloc(n * sizeof(double*));
	if (ans == NULL)
		return NULL;
	for (i = 0; i < n; i++) {
		ans[i] = (double*) malloc(n * sizeof(double));
		if (ans[i] == NULL)
			return NULL;
	}
	for (x = 0; x < n; x++)
		for (y = 0; y < n; y++)
			ans[x][y] = 0;
	for (x = 0; x < n; x++) {
		ans[x][x] = 1;
	}
	//确定最大层数和最小层数
	double rand_max = log(n) / log(2);
	double rand_min = (double) (((int) (log2(l)))) + 1;
	//声明循环中的变量
	double nn = 0;
	double *p1_0 = NULL, *p2_0 = NULL;
	double **p1 = NULL, **p2 = NULL;
	double **w1 = NULL;
	//循环开始
	for (i = rand_min; i <= rand_max; i++) {
		nn = pow(2, i);
		//为循环中用到的变量声明内存空间，并清零
		p1_0 = (double*) malloc(nn * sizeof(double));
		p2_0 = (double*) malloc(nn * sizeof(double));
		memset(p1_0, 0, sizeof(double));
		memset(p2_0, 0, sizeof(double));
		//构造向量
		for (x = 0; x < 16; x++)
			p1_0[x] = h[x];
		for (x = 0; x < 16; x++)
			p2_0[x] = g[x];
		//为向量圆周移位准备内存空间，如果分配失败返回NULL
		p1 = (double**) malloc((nn / 2) * sizeof(double*));
		if (p1 == NULL)
			return NULL;
		p2 = (double**) malloc((nn / 2) * sizeof(double*));
		if (p2 == NULL)
			return NULL;
		for (x = 0; x < nn / 2; x++) {
			p1[x] = (double*) malloc(nn * sizeof(double));
			if (p1[x] == NULL)
				return NULL;
		}
		for (x = 0; x < nn / 2; x++) {
			p2[x] = (double*) malloc(nn * sizeof(double));
			if (p2[x] == NULL)
				return NULL;
		}
		for(i=0;i<nn/2;i++)
					for(j=0;j<nn;j++){
						p1[i][j]=0;
						p2[i][j]=0;
					}
		//开始向量圆周移位
		for (j = 0; j < nn / 2; j++) {
			for (x = 0; x < nn / 2; x++) {
				p1[j][x] = p1_0[x];
				p2[j][x] = p2_0[x];
			}
			temp1 = p1_0[(int) nn - 1];
			temp2 = p1_0[(int) nn - 2];
			temp3 = p2_0[(int) nn - 1];
			temp4 = p2_0[(int) nn - 2];
			for (y = 31; y >= 2; y--) {
				p1_0[y] = p1_0[y - 2];
				p2_0[y] = p2_0[y - 2];
			}
			p1_0[0] = temp2;
			p1_0[1] = temp1;
			p2_0[0] = temp4;
			p2_0[1] = temp3;
		}
		free(p1_0);
		free(p2_0);

		//构造正交矩阵
		//为临时用矩阵申请空间,并赋值
		w1 = (double**) malloc(nn * sizeof(double*));
		if (!w1)
			return NULL;
		for (x = 0; x < nn; x++) {
			w1[x] = (double*) malloc(nn * sizeof(double));
			if (w1[x] == NULL)
				return NULL;
		}
		for (x = 0; x < nn / 2; x++) {
			for (y = 0; y < nn; y++) {
				w1[x][y] = p1[x][y];
			}
		}
		for (x = nn / 2; x < nn; x++) {
			for (y = 0; y < nn; y++) {
				w1[x][y] = p2[x - (int) (nn / 2)][y];
			}
		}
		ans = matrixmultiplic(w1, n, ans, n,n);
		//clear p1;clear p2;
		clear(p1, nn / 2);
		clear(p2, nn / 2);
	}
	return ans;
}

