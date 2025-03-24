// MyCpp.h: interface for the CMyCpp class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MYCPP_H__C63EB901_D617_4F46_96FE_5B8B28A71665__INCLUDED_)
#define AFX_MYCPP_H__C63EB901_D617_4F46_96FE_5B8B28A71665__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <stdio.h>


//////////////////////////////////////////////////////////////////////////
//       �Գƾ����±���㺯��

#define ik (i>=k)? i*(i+1)/2+k : k*(k+1)/2+i
#define kj (k>=j)? k*(k+1)/2+j : j*(j+1)/2+k




//  Ȩ�������㺯��
double Calculate_q(double *Q,double *F,int t);
//  Ȩ�������㺯��
double Calculate_q(double Q[],double F[],int Fin[],int n);




void MyBreak(const char* fmt, ...);

int ij(int i,int j);



bool inverse(double a[],int n);

/*
bool inverse3(double A[],int n);
bool inverse4(double A[],int n);
*/


// Ȩ���󴫲�����
void Calculate_BQBT(double B[],double Q[],int r,int n,double N[]);



//  ���ȷ�����д��(double��)�ǶȻ�Ϊ����ֵ
double dms_rad(double a);

//  ���ǶȵĻ���ֵ��Ϊ�ȷ�����д�ĽǶȣ�double �ͣ� 
double rad_dms(double a);


void PrintM(FILE *fp,double A[],int size, int t,char* fmt,
			char* title=NULL,bool IsLabel=true);

void  PrintM2(FILE* fp, double M[], int n,int t,const char *fmt,
			  const char* title=NULL,bool IsLabel=true);

void  PrintEquation(FILE* fp, double A[], double b[],
					int n, int m,char *fmt,char* title=NULL);


#define IGG1  0
#define IGG3  1
#define Huber 2

double Wi(int fname, double v, double k0, double k1);

double Median(double pp[],int n,bool IsAbs);//��λ������

void calculateTransformationParameters(const int n, const double* srcPoints, const double* dstPoints, double& dx, double& dy, double& scale, double& rotation);//// �����Ĳ���ת�� 
void transformPoint(const double x, const double y, double& rex, double& rey, double dx, double dy, double scale, double rotation);// ʹ���Ĳ�����������ת�� 


#endif // !defined(AFX_MYCPP_H__C63EB901_D617_4F46_96FE_5B8B28A71665__INCLUDED_)
