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
double *OMP(double* ss, double** rr, int height, int width, int N)
{

    double* hat_y=NULL;
    double** r=(double*) malloc(height * sizeof(double));
    double** s=(double**)malloc(height*sizeof(double*));
    double** Aug_t =(double**)malloc(height*sizeof(double*));
    double* r_n = (double*) malloc(height * sizeof(double));
    double* product = (double*) malloc(N * sizeof(double));
    double** aug_y = NULL,**Aug_tt=NULL;
    double** temp = NULL, **temp1 = NULL, **temp3 = NULL, **temp4 = NULL;
    int pos_array[256]= {0};
    int M = height, i,j;

    for(i=0; i<height; i++)
    {
        s[i]=(double*)malloc(sizeof(double));
        r[i]=(double*)malloc(width*sizeof(double));
        Aug_t[i]=(double*)malloc(width*sizeof(double));
        s[i][0]=ss[i];
    }
    for(i=0; i<height; i++)
        for(j=0; j<width; j++)
            r[i][j]=rr[i][j];
    int times, col, pos;

//	int count=0;//Aug_t的长度为times
    for (i = 0; i < height; i++)
        r_n[i] = s[i][0];

    for (times = 0; times < M; times++)
    {
        printf("%d start OK\n", times+1);
        for (col = 0; col < N; col++)
        {
//////    for col=1:N;                                  %  恢复矩阵的所有列向量
//////        product(col)=abs(T(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值)
//////    end
            double* colu=getcolumn(r, height, col);
            product[col] = fabs(innerproduct(colu, r_n, height));
        }
        printf("%d product(col)=abs(T(:,col)'*r_n); OK\n", times+1);
//////    [val,pos]=max(product);                       %  最大投影系数对应的位置
//        testoutnn(product,256,"product.txt");
        pos = max(product, N);
        printf("%d [val,pos]=max(product); %d OK\n", times+1,pos);
//////    Aug_t=[Aug_t,T(:,pos)];                       %  矩阵扩充
        expandmatrix(Aug_t, getcolumn(r, height, pos), height, times);
        printf("%d Aug_t=[Aug_t,T(:,pos)]; OK\n", times+1);
//        testoutn(Aug_t,height,times+1,"Aug_t.txt");
//////    T(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零）
        formatcolumn(r, height, pos, 0);
        printf("%d T(:,pos)=zeros(M,1); OK\n", times+1);
//        testoutn(r,height,width,"r.txt");
//////    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  最小二乘,使残差最小
//////    r_n=s-Aug_t*aug_y;                            %  残差
//////    pos_array(times)=pos;                         %  纪录最大投影系数的位置
//////
//////    if (abs(aug_y(end))^2/norm(aug_y)<0.05)       %  自适应截断误差（***需要调整经验值）
//////        break;
//////    end

        //aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;
        Aug_tt=matrixturn(Aug_t, height, times+1);
//        printf("%d matrixturn OK\n", times+1);
        temp = matrixmultiplic(Aug_t, times+1,Aug_tt,times+1,M);
//        printf("%d matrixmultiplic OK\n", times+1);
        temp1 = inversematrix(temp, times+1);
//        printf("%d inversematrix OK\n", times+1);
        temp3 = matrixmultiplic(Aug_tt, height, temp1, times+1,times+1);
//        printf("%d matrixmultiplic OK\n", times+1);
        aug_y = matrixmultiplic(s, 1, temp3, times+1,M);	//aug_y好多行一列
//        testoutn(aug_y,times+1,1,"aug_y.txt");
        printf("%d aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s; OK\n", times+1);
        //r_n=s-Aug_t*aug_y;
        temp4 = matrixmultiplic(aug_y, 1, Aug_t, height,times+1);
//        testoutn(temp4,height,1,"temp4.txt");
//        printf("%d matrixmultiplic OK\n", times+1);
        r_n = vectorsub(s, temp4, height);
        printf("%d r_n=s-Aug_t*aug_y; OK\n", times+1);
        //pos_array(times)=pos;
        vectoradd(pos_array, times, pos);
        printf("%d pos_array(times)=pos; OK\n\n", times+1);
        //if (abs(aug_y(end))^2/norm(aug_y)<0.05)
        if(pow(fabs(aug_y[times][0]),2)/norm(aug_y,times)<0.05) break;
        clear(temp,times+1);
    }
    hat_y=finish(pos_array,aug_y,times,N);
    for(i=0; i<height; i++)
        free(r[i]);
    free(r);
    return hat_y;
}
