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
#include "DWT.h"
#include "matrixop.h"
#include "OMPOP.h"

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

void testout(double** x, int height, int width, char*str);	//将指定矩阵输出至选中的文件中
void writebmp(FILE* fo, double** pixel, int height, int width);
void readbmp(FILE* fp, double** pixel, int height, int width);	//读取图像
void writebmp(FILE* fo, double** pixel, int height, int width);
double** createsimplingmatrix(int wid, int height);	//生成观测矩阵
double error(double** piexl, double** ans, int height, int width);
//double** observation(double**, int, double**, int);	//进行观测
double gaussrand();	//生成高斯矩阵

int main() {
	FILE *fp,  *fo;
	BMPFILEHEADER fileHeader;
	BMPINF infoHeader;
	long width, height;
	int i;
	int temp[256*4];
//		long offset,bmpImageSize, bytesPerPixel, size, bitCount;
//		WORD c;
	double **pixel = NULL;
	double **simplingmatrix = NULL;
	double **compressed_matrix = NULL;
	double **ww = NULL, **wwturn = NULL;
	double** temp1=NULL;
	double **ans = NULL;

	if ((fp = fopen("lena256.bmp", "rb")) == NULL) {
		printf("Cann't open the file!\n");
		exit(0);
	}
	if ((fo = fopen("OK.bmp", "wb")) == NULL) {
		printf("Cann't open the file!\n");
		exit(0);
	}
	fseek(fp, 0, 0);
	fread(&fileHeader, sizeof(fileHeader), 1, fp);
	fwrite(&fileHeader, sizeof(fileHeader), 1, fo);
	fread(&infoHeader, sizeof(infoHeader), 1, fp);
	fwrite(&infoHeader, sizeof(infoHeader), 1, fo);
	fread(&temp,sizeof(int),256*4,fp);
	fwrite(&temp,sizeof(int),256*4,fo);
	//fseek(fp, 256 * 4, SEEK_CUR);

	//计算并输出位图数据的偏移量，图像的大小，宽度和高度，每个像素点所占的字节
	width = infoHeader.bWidth;
	height = infoHeader.bHeight;
//		size = fileHeader.bSize;
//		offset = fileHeader.bOffset;
//		bmpImageSize = infoHeader.bmpImageSize;
//		bitCount = infoHeader.bBitCount;
//		bytesPerPixel = infoHeader.bBitCount / 8;
	pixel = creatmatrix(height, width);
	if (pixel == NULL)
		return 0;
	readbmp(fp, pixel, height, width);
	simplingmatrix = createsimplingmatrix(190, height);
	compressed_matrix = matrixmultiplic(pixel, width, simplingmatrix, 190,width);
	ww = creatmatrix(height, height);
	ww = DWT(height);
	wwturn = matrixturn(ww, height, height);
	compressed_matrix = matrixmultiplic(wwturn, width, compressed_matrix,
			190,width);
	simplingmatrix = matrixmultiplic(wwturn, width, simplingmatrix, 190,width);
	ans = creatmatrix(height, width);
	/*for (i = 0; i < width; i++) {
		double* rec = OMP(getcolumn(compressed_matrix, 190, i),
				simplingmatrix, 190, width, height);
		setcolum(ans, i, rec, height);
	}
	//小波反变换
	temp1=matrixmultiplic(matrixturn(ans,height,width),width,ww,height,width);
	ans=matrixmultiplic(ans,width,temp1,height,width);*/
	//输出结果
	writebmp(fo,ans,height,width);
	//误差

	//测试用部分
	//testout(pixel,height,width,"test.txt");
	//testout(compressed_matrix,190,width,"compressed_matrix.txt");
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
void readbmp(FILE* fp, double** pixel, int height, int width) {
	int i, j;
	for (i = height - 1; i >= 0; i--) {
		for (j = 0; j < width; j++) {
			pixel[i][j] = (double) fgetc(fp);
		}
	}
	fclose(fp);
}
double** createsimplingmatrix(int wid, int height) {
	double **phi = NULL;
	int i, j;
	phi = (double**) malloc(height * sizeof(double*));
	if (phi == NULL)
		return NULL;
	for (i = 0; i < height; i++) {
		phi[i] = (double*) malloc(wid * sizeof(double));
		if (phi[i] == NULL)
			return NULL;
	}
	for (i = 0; i < height; i++) {
		for (j = 0; j < wid; j++) {
			phi[i][j] = gaussrand();
		}
	}
	return phi;
}
double gaussrand() {
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if (phase == 0) {
		do {

			double U1 = (double) rand() / RAND_MAX;
			double U2 = (double) rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);

		X = (double) (V1 * sqrt(-2 * log(S) / S));
	} else
		X = (double) (V2 * sqrt(-2 * log(S) / S));

	phase = 1 - phase;

	return X;
}
void writebmp(FILE* fo, double** pixel, int height, int width) {
	int i, j;
	int temp[height][width];
	for(i=0;i<height;i++)
		for(j=0;j<width;j++)
			temp[i][j]=(int)pixel[i][j];
		for (i = height - 1; i >= 0; i--) {
			for (j = 0; j < width; j++) {
				fwrite(&temp[i][j],sizeof(int),1,fo);
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
/*double** observation(double** pixel, int width, double** simplingmatrix,
 int sheight) {
 double** temp = NULL;
 int i, j, k;
 //	动态分配空间，为观测后的矩阵准备存储空间
 temp = (double**) malloc(sheight * sizeof(double*));
 if (temp == NULL)
 return NULL;
 for (i = 0; i < sheight; i++) {
 temp[i] = (double*) malloc(width * sizeof(double));
 if (temp[i] == NULL)
 return NULL;
 }
 //	矩阵乘法，进行观测
 for (i = 0; i < sheight; i++) {
 for (j = 0; j < width; j++) {
 temp[i][j] = 0;
 for (k = 0; k < width; k++) {
 temp[i][j] += simplingmatrix[i][k] * pixel[k][j];
 }
 }
 }
 return temp;
 }*/
