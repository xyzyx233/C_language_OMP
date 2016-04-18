/*
 ============================================================================
 Name        : ���м���_����.c
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

/*1������ת�� ֱ����
 2������˷� ����Ϊ��������
 3��ָ����С������� ����Ϊ����ĳ���
 4���ƶ���С��ȫ��1�ľ��� ����Ϊ����ĳ���
 5���������ĳ��� ����Ϊ���󣬷�����Ϊ1*2�����飬��һ���������������ڶ�������������
 ��Щ������λͼ��ÿ��Ԫ��Ϊһ�����ص�*/
/*
 ����λͼ�������ļ��е�����˳���Ǵ����½ǵ����Ͻǣ�����Ϊ�������еġ�
 */

/*����WORDΪ�����ֽڵ�����*/
typedef unsigned short WORD;
/*����DWORDΪe�ĸ��ֽڵ�����*/
typedef unsigned long DWORD;
/*λͼ�ļ�ͷ*/
typedef struct BMP_FILE_HEADER {
	WORD bType; /* �ļ���ʶ�� */
	DWORD bSize; /* �ļ��Ĵ�С */
	WORD bReserved1; /* ����ֵ,��������Ϊ0 */
	WORD bReserved2; /* ����ֵ,��������Ϊ0 */
	DWORD bOffset; /* �ļ�ͷ�����ͼ������λ��ʼ��ƫ���� */
} BMPFILEHEADER;
/*λͼ��Ϣͷ*/
typedef struct BMP_INFO {
	DWORD bInfoSize; /* ��Ϣͷ�Ĵ�С */
	DWORD bWidth; /* ͼ��Ŀ�� */
	DWORD bHeight; /* ͼ��ĸ߶� */
	WORD bPlanes; /* ͼ���λ���� */
	WORD bBitCount; /* ÿ�����ص�λ�� */
	DWORD bCompression; /* ѹ������ */
	DWORD bmpImageSize; /* ͼ��Ĵ�С,���ֽ�Ϊ��λ */
	DWORD bXPelsPerMeter; /* ˮƽ�ֱ��� */
	DWORD bYPelsPerMeter; /* ��ֱ�ֱ��� */
	DWORD bClrUsed; /* ʹ�õ�ɫ���� */
	DWORD bClrImportant; /* ��Ҫ����ɫ�� */
} BMPINF;
/*��ɫ��*/
typedef struct RGB_QUAD {
	WORD rgbBlue; /* ��ɫǿ�� */
	WORD rgbGreen; /* ��ɫǿ�� */
	WORD rgbRed; /* ��ɫǿ�� */
	WORD rgbReversed; /* ����ֵ */
} RGBQUAD;

void testout(double** x, int height, int width, char*str);	//��ָ�����������ѡ�е��ļ���
void writebmp(FILE* fo, double** pixel, int height, int width);
void readbmp(FILE* fp, double** pixel, int height, int width);	//��ȡͼ��
void writebmp(FILE* fo, double** pixel, int height, int width);
double** createsimplingmatrix(int wid, int height);	//���ɹ۲����
double error(double** piexl, double** ans, int height, int width);
//double** observation(double**, int, double**, int);	//���й۲�
double gaussrand();	//���ɸ�˹����

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

	//���㲢���λͼ���ݵ�ƫ������ͼ��Ĵ�С����Ⱥ͸߶ȣ�ÿ�����ص���ռ���ֽ�
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
	//С�����任
	temp1=matrixmultiplic(matrixturn(ans,height,width),width,ww,height,width);
	ans=matrixmultiplic(ans,width,temp1,height,width);*/
	//������
	writebmp(fo,ans,height,width);
	//���

	//�����ò���
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
 //	��̬����ռ䣬Ϊ�۲��ľ���׼���洢�ռ�
 temp = (double**) malloc(sheight * sizeof(double*));
 if (temp == NULL)
 return NULL;
 for (i = 0; i < sheight; i++) {
 temp[i] = (double*) malloc(width * sizeof(double));
 if (temp[i] == NULL)
 return NULL;
 }
 //	����˷������й۲�
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
