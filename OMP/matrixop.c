/*
 * matrixop.c
 *
 *  Created on: 2016年4月5日
 *      Author: asus
 */

#include <math.h>
#include  <stdio.h>
#include <stdlib.h>

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
double** matrixmultiplic(double** matrix1, int width, double** matrix2,
		int height, int wid2) {
	int i, j, k;
	double** ans = NULL;
	//double** ans = creatmatrix(height,width); 这种写法会报错 conflicting types for 在声明前使用
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
			for (k = 0; k < wid2; k++) {
				ans[i][j] += matrix2[i][k] * matrix1[k][j];
			}
		}
	}
	return ans;
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
double* createvector(int n) {
	double* ans = (double*) malloc(n * sizeof(double));
	int i = 0;
	if (!ans)
		return NULL;
	for (i = 0; i < n; i++)
		ans[i] = 0;
	return ans;
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
double* getrow(double** matrix, int column, int row) {
	int i;
	double* vactor = createvector(column);
	for (i = 0; i < row; i++) {
		matrix[i][column]=vactor[i];
	}
	return vactor;
}
void setcolum(double** matrix, int column, double* vactor, int row) {
	int i;
	for (i = 0; i < row; i++) {
		matrix[column][i] = vactor[i];
	}
}
double innerproduct(double* a, double* b, int len) {
	int i;
	double sum = 0;
	for (i = 0; i < len; i++) {
		sum += a[i] * b[i];
	}
	return sum;
}

int maxx(double* vector, int num) {
	int imax=0,i=0;
	for(i=0;i<num;i++)
	        if(vector[imax]<vector[i])
	            imax=i;
	return imax;
}
void formatcolumn(double** matrix, int height, int column, int val) {
	int i;
	for (i = 0; i < height; i++) {
		matrix[i][column] = val;
	}
}
void expandmatrix(double **ma1, double *ve, int height, int width) {
	int i;
	for (i = 0; i < height; i++) {
		ma1[i][width] = ve[i];
	}
}
void inversematrix(double **A_, int n,double** out) { //计算每一行每一列的每个元素所对应的余子式，组成A*

    double **A = (double**)malloc(n * sizeof(double*));
    double **L = (double**)malloc(n * sizeof(double*));
    double **U = (double**)malloc(n * sizeof(double*));
    double **u = (double**)malloc(n * sizeof(double*));
    double **r = (double**)malloc(n * sizeof(double*));
    //double **out;
    int i,j,k,x;
    double s;
    //out = (double**)malloc(n*sizeof(double*));
    for( i = 0; i < n; i++)
    {
        A[i] = (double*)malloc(n*sizeof(double));
        for( j = 0; j < n; j++)
        {
            A[i][j] = A_[i][j];
        }
    }
    for( i = 0; i < n; i++)
    {
        L[i] = (double*)malloc((i+1)*sizeof(double));
        U[i] = (double*)malloc(n*sizeof(double));
        u[i] = (double*)malloc(n*sizeof(double));
        r[i] = (double*)malloc(n*sizeof(double));
        //out[i] = (double*)malloc(n*sizeof(double));
        L[i][i] = 1;
    }
	for(i=0;i<n;i++)
		for (j = 0;j < n;j++) {
			u[i][j] = 0;
			r[i][j] = 0;
		}
    //直接三角分解
    //首先计算矩阵L的第1列，矩阵U的第1行
    for( i = 0; i < n; i++)
    {
        U[0][i] = A[0][i];
        if(i != 0)
        {
            L[i][0] = A[i][0]/U[0][0];
        }
    }
    //递推求矩阵L第r行和矩阵U第r列
    for(x = 1; x < n; x++)
    {
        for(int i = x; i < n; i++)
        {
            double tmp1 = 0;
            double tmp2 = 0;
            for(int k = 0; k < x; k++)
            {
                tmp1 += L[x][k]*U[k][i];
                if(i > x)
                    tmp2 += L[i][k]*U[k][x];
            }
            U[x][i] = A[x][i] - tmp1;
            if(i > x)
                L[i][x] = (A[i][x] - tmp2)/U[x][x];
        }
    }
/////////////////////求L和U矩阵的逆//////////////////////////////////////////
        for (i=0; i<n; i++) /*求矩阵U的逆 */
        {
            u[i][i]=1/U[i][i];//对角元素的值，直接取倒数
            for (k=i-1; k>=0; k--)
            {
                s=0;
                for (j=k+1; j<=i; j++)
                    s=s+U[k][j]*u[j][i];
                u[k][i]=-s/U[k][k];//迭代计算，按列倒序依次得到每一个值，
            }
        }
        for (i=0; i<n; i++) //求矩阵L的逆
        {
            r[i][i]=1; //对角元素的值，直接取倒数，这里为1
            for (k=i+1; k<n; k++)
            {
                for (j=i; j<=k-1; j++)
                    r[k][i]=r[k][i]-L[k][j]*r[j][i];   //迭代计算，按列顺序依次得到每一个值
            }
        }
//////////将r和u相乘，得到逆矩阵
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                out[i][j]=0;
            }
        }
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                for(k=0; k<n; k++)
                {
                    out[i][j]+=u[i][k]*r[k][j];
                }
            }
        }
        for(i=0;i<n;i++){
            free(A[i]);
            free(U[i]);
            free(u[i]);
            free(r[i]);
            free(L[i]);
        }
        free(A);
        free(L);
        free(U);
        free(u);
        free(r);
		A = NULL;
		L = NULL;
		U = NULL;
		u = NULL;
		r = NULL;

        //return out;
}
double* vectorsub(double** vec1, double** vec2, int num) {
	double* ans = (double*) malloc(num * sizeof(double));
	int i;
	for (i = 0; i < num; i++)
		ans[i] = vec1[i][0] - vec2[i][0];
	return ans;
}
void vectoradd(int* vec, int num, int n) {
	vec[num] = n;
}
double norm(double* vec,int num){
	double ans=0;
	int i;
	for(i=0;i<num;i++)
		ans+=vec[i]*vec[i];
	return sqrt(ans);
}
double* finish(int* pos, double** aug, int num,int n){
	double* ans=createvector(n);
	int i;
	for(i=0;i<num;i++){
		ans[pos[i]]=aug[i][0];
	}
	return ans;
}
void Reverse(double *arr, int b, int e)
{
    for(; b < e; b++, e--)
    {
        double temp = arr[e];
        arr[e] = arr[b];
        arr[b] = temp;
    }
}

void RightShift(double *arr, int N, int k)
{
    k %= N;
    Reverse(arr, 0, N-k-1);
    Reverse(arr, N-k, N-1);
    Reverse(arr, 0, N-1);
}
