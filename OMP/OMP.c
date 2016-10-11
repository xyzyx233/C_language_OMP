/*
 ============================================================================
 Name        : 并行计算_串行.c
 Author      : eric
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <memory.h>
#include "matrixop.h"
#include "OMPOP.h"

#define RR 190
#define M 256
#define N 256

#pragma pack(2)

/*1、矩阵转置 直接用
 2、矩阵乘法 参数为两个矩阵
 3、指定大小的零矩阵 参数为矩阵的长宽
 4、制定大小的全是1的矩阵 参数为矩阵的长宽
 5、求出矩阵的长宽 参数为矩阵，返回量为1*2的数组，第一个用来存行数，第二个用来存列数
 这些矩阵都是位图，每个元素为一个像素点*/
/*
 所以位图数据在文件中的排列顺序是从左下角到右上角，以行为主序排列的。
 */

/*定义WORD为两个字节的类型*/
typedef unsigned short WORD;
/*定义DWORD为e四个字节的类型*/
typedef unsigned long DWORD;
/*位图文件头*/
typedef struct BMP_FILE_HEADER {
	WORD bType; /* 文件标识符 */
	DWORD bSize; /* 文件的大小 */
	WORD bReserved1; /* 保留值,必须设置为0 */
	WORD bReserved2; /* 保留值,必须设置为0 */
	DWORD bOffset; /* 文件头的最后到图像数据位开始的偏移量 */
} BMPFILEHEADER;
/*位图信息头*/
typedef struct BMP_INFO {
	DWORD bInfoSize; /* 信息头的大小 */
	DWORD bWidth; /* 图像的宽度 */
	DWORD bHeight; /* 图像的高度 */
	WORD bPlanes; /* 图像的位面数 */
	WORD bBitCount; /* 每个像素的位数 */
	DWORD bCompression; /* 压缩类型 */
	DWORD bmpImageSize; /* 图像的大小,以字节为单位 */
	DWORD bXPelsPerMeter; /* 水平分辨率 */
	DWORD bYPelsPerMeter; /* 垂直分辨率 */
	DWORD bClrUsed; /* 使用的色彩数 */
	DWORD bClrImportant; /* 重要的颜色数 */
} BMPINF;
/*彩色表*/
typedef struct RGB_QUAD {
	WORD rgbBlue; /* 蓝色强度 */
	WORD rgbGreen; /* 绿色强度 */
	WORD rgbRed; /* 红色强度 */
	WORD rgbReversed; /* 保留值 */
} RGBQUAD;

double anss[M][256];
unsigned char ans1[M][256];
double ww1[M][256];
double wwturn1[M][256];
double tempx[M][256];

void testout(double** x, int height, int width, char*str);	//将指定矩阵输出至选中的文件中
void writebmp(FILE* fo, double** pixel, int height, int width);
void readbmp(FILE* fp, double** pixel);	//读取图像
void writebmp(FILE* fo, unsigned int pixel[][N], int height, int width);
double** createsimplingmatrix(int wid, int height);	//生成观测矩阵
double error(double** piexl, double** ans, int height, int width);
//double** observation(double**, int, double**, int);	//进行观测
void matrixmultiplicxx(double matrix1[][N], int width, double matrix2[][256],
	int height, int wid2, double ans[][N]) {
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
void testoutcol(double* x, int height, char*str) {
	int i;
	FILE *fto = fopen(str, "wb");
	for (i = 0; i < height; i++) {
			fprintf(fto, "%f\t", x[i]);
	}
	fclose(fto);
}
double gaussrand();	//生成高斯矩阵

int main() {
	FILE *fp,  *fo, *fp1;
	FILE *fp2, *fp3;
	BMPFILEHEADER fileHeader;
	BMPINF infoHeader;
	long width, height;
	int i,j;
	int temp[256 *4];
//		long offset,bmpImageSize, bytesPerPixel, size, bitCount;
//		WORD c;
	double **pixel=NULL;
	double **simplingmatrix=NULL;
	double **compressed_matrix = NULL;
	double **ww =NULL, **wwturn=NULL;
	double **ans = NULL;
	double **temp1;

	if ((fp = fopen("lena256.bmp", "rb")) == NULL) {
		printf("Cann't open lena256!\n");
		exit(0);
	}
	if ((fp1 = fopen("DWT.txt", "rb")) == NULL) {
		printf("Cann't open DWT.txt!\n");
		exit(0);
	}
	if ((fo = fopen("OK.bmp", "wb")) == NULL) {
		printf("Cann't open OK.bmp!\n");
		exit(0);
	}
	fseek(fp, 0, 0);
	fread(&fileHeader, sizeof(fileHeader), 1, fp);
	fwrite(&fileHeader, sizeof(fileHeader), 1, fo);
	fread(&infoHeader, sizeof(infoHeader), 1, fp);
	fwrite(&infoHeader, sizeof(infoHeader), 1, fo);
	fread(&temp,sizeof(unsigned int),256,fp);
	fwrite(&temp,sizeof(unsigned int),256,fo);

	//计算并输出位图数据的偏移量，图像的大小，宽度和高度，每个像素点所占的字节
	width = infoHeader.bWidth;
	height = infoHeader.bHeight;
	pixel = creatmatrix(height, width);
	if (pixel == NULL)
		return 0;
	readbmp(fp, pixel);

	simplingmatrix = createsimplingmatrix(RR, width);
	compressed_matrix = matrixmultiplic(pixel, width, simplingmatrix, RR,width);

	ww = creatmatrix(M, N);
//	ww = DWT(height);
    for(i=0;i<M;i++)
        for(j=0;j<N;j++)
        fscanf(fp1,"%lf",&ww[i][j]);

	wwturn = matrixturn(ww, height, height);
	for(i=0;i<M;i++)
		for (j = 0;j < N;j++) {
			ww1[i][j] = ww[i][j];
			wwturn1[i][j] = wwturn[i][j];
		}
	temp1 = matrixmultiplic(wwturn, width, compressed_matrix,
		RR,width);
	compressed_matrix=temp1;

	temp1 = matrixmultiplic(wwturn, width, simplingmatrix, RR,width);
	simplingmatrix=temp1;
    ans = creatmatrix(height, width);

//testout(compressed_matrix, 190, width, "compressed_matrix.txt");
//testout(simplingmatrix, 190, width, "simplingmatrix.txt");
	//simplingmatrix = creatmatrix(190, 256);
	//compressed_matrix = creatmatrix(190, 256);
	//if ((fp2 = fopen("s.txt", "rb")) == NULL) {
	//	printf("Cann't open DWT.txt!\n");
	//	exit(0);
	//}
	//for (i = 0;i<190;i++)
	//	for (j = 0;j<256;j++)
	//		fscanf(fp2, "%lf", &simplingmatrix[i][j]);
	//if ((fp3 = fopen("c.txt", "rb")) == NULL) {
	//	printf("Cann't open DWT.txt!\n");
	//	exit(0);
	//}
	//for (i = 0;i<190;i++)
	//	for (j = 0;j<256;j++)
	//		fscanf(fp3, "%lf", &compressed_matrix[i][j]);
	printf("OMP start\n");
    for (i = 0; i < width; i++)
    {
        double* col=getcolumn(compressed_matrix, RR, i);
        double* rec =OMP(col, simplingmatrix, RR, width, height);
        setcolum(ans, i, rec, height);
		printf("%d OMP OK\n",i+1);
    }
	for (i = 0;i < M;i++)
		for (j = 0;j < N;j++)
			anss[i][j] = ans[i][j];
	//小波反变换
	matrixmultiplicxx(anss,width,wwturn1,height,width,tempx);
	matrixmultiplicxx(ww1,width,tempx,height,width,anss);
	printf("小波变换完成\n");
	for (i = 0;i < M;i++)
		for (j = 0;j < N;j++)
			ans[i][j] = anss[i][j];
	testout(ans, height, width, "ans.txt");
	printf("ans 输出完成\n");

	for (i = 0;i < M;i++) {
		for (j = 0;j < N;j++) {
			if (anss[i][j] > 255) {
				ans1[i][j] = 255;
			}
			else if (anss[i][j] < 0) {
				ans1[i][j] = 0;
			}
			else {
				ans1[i][j] = anss[i][j];
			}
		}
	}
	printf("准备输出图像\n");

	//FILE *fto = fopen("ans1.txt", "wb");
	//for (i = 0; i < height; i++) {
	//	for (j = 0; j < width; j++) {
	//		fprintf(fto, "%d\t", ans1[i][j]);
	//	}
	//	fprintf(fto, "\n");
	//}
	//fclose(fto);
//	//输出结果
	//for (i = 0;i < 256;i++)
	//	for (j = 0;j < 256;j++)
	//		ans1[i][j] = pixel[i][j];
	writebmp(fo,ans1,height,width);

	//测试用部分
	//testout(ans, height, width, "ans.txt");
	//testout(pixel,height,width,"pixel.txt");
	return 0;
}
void testout(double** x, int height, int width, char*str) {
	int i, j;
	FILE *fto = fopen(str, "wb");
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			fprintf(fto, "%f\t", x[i][j]);
		}
		fprintf(fto, "\n");
	}
	fclose(fto);
}
void readbmp(FILE* fp, double** pixel) {
	int i, j;
	for (i = 256 - 1; i >= 0; i--) {
		for (j = 0; j < 256; j++) {
			pixel[i][j] = (double) fgetc(fp);
		}
	}
	fclose(fp);
}
double** createsimplingmatrix(int height, int wid) {
	double **phi = NULL;
	int i, j;
	double x, p[190][256];
	phi = (double**) malloc(height * sizeof(double*));
	if (phi == NULL)
		return NULL;
	for (i = 0; i < height; i++) {
		phi[i] = (double*) malloc(wid * sizeof(double));
		if (phi[i] == NULL)
			return NULL;
	}
	for(i=0;i<height;i++)
		for(j=0;j<wid;j++){
			phi[i][j]=0;
			p[i][j]=0;
		}
	srand((int)time(NULL));
	for (i = 0; i < height; i++) {
		for (j = 0; j < wid; j++) {
//			while(1){
				x=gaussrand();
//				if(x>=-10||x<=10){
//					break;
//				}
//			}
			p[i][j] = x;
		}
	}
	for(i=0;i<height;i++)
		for(j=0;j<wid;j++){
			phi[i][j]=p[i][j];
		}
	return phi;
}
double gaussrand() {
    static double U, V;
    static int phase = 0;
    double z;

    if(phase == 0)
    {
         U = (rand() + 1.1) / (RAND_MAX +2.);
         V = rand() / (RAND_MAX + 1.);
         z = sqrt(-1 * log(U))* sin(2 * 3.141592654 * V);
    }
    else
    {
         z = sqrt(-2 * log(U)) * cos(2 * 3.141592654 * V);
    }

    phase = 1 - phase;
    return z;
}
void writebmp(FILE* fo, unsigned char pixel[][N], int height, int width) {
	int i, j;
		for (i = height - 1; i >= 0; i--) {
			for (j = 0; j < width; j++) {
				fwrite(&pixel[j][i],sizeof(unsigned char),1,fo);
			}
		}
		fclose(fo);
}
double error(double** piexl, double** ans, int height, int width){
	double sum=0, psnr;
	int i,j;
	for(i=0;i<height;i++)
		for(j=0;j<width;j++)
			sum+=pow(fabs(ans[i][j]-piexl[i][j]),2);
	psnr=10*log10(255*255/(sum/height/width));
	return psnr;
}
