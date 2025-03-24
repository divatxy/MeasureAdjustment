// MyCpp.cpp: implementation of the CMyCpp class.
//
//////////////////////////////////////////////////////////////////////
#define _CRT_SECURE_NO_WARNINGS

//#include "stdafx.h"
#include "stdio.h"
#include <iostream>
#include <stdarg.h>
#include <stdlib.h>
#include "math.h"
#include "public.h"
#include <vector> 
#include <cmath> 
#include <Eigen/Dense> // 需要Eigen库 
#define PI 3.14159265358979323846

//////////////////////////////////////////////////////////////////////////
//   显示提示信息
void MyBreak(const char* fmt, ...)
{
	char buffer[256];
	va_list argptr;
	va_start(argptr, fmt);
	vsprintf(buffer, fmt, argptr);
	va_end(argptr);
	
#ifdef VC_EXTRALEAN
	AfxMessageBox(buffer);
#else
	printf(buffer);
	getchar();
#endif // VC_EXTRALEAN
	
}



//////////////////////////////////////////////////////////////////////////
//       对称矩阵下标计算函数
int ij(int i,int j)
{
	return (i>=j)? i*(i+1)/2+j :j*(j+1)/2+i;
}


//////////////////////////////////////////////////////////////////////////
//  权倒数计算函数
double Calculate_q(double *Q,double *F,int t)
{
	double q=0.0;
	for(int i=0;i<t;i++)
		for(int j=0;j<t;j++) 
			q+=Q[ij(i,j)]*F[i]*F[j];
		
	return q;
		
}



//////////////////////////////////////////////////////////////////////////
//  权倒数计算函数(系数系列进行了压缩)
double Calculate_q(double Q[],double F[],int Fin[],int n)
{
	double q=0.0;
	for(int k=0;k<n;k++)
	{
		int i=Fin[k];
		for(int s=0;s<n;s++) 
		{
			int j=Fin[s];
			q+=Q[ij(i,j)]*F[k]*F[s];
		}
	}
	return q;
	
}




//////////////////////////////////////////////////////////////////////////
//  将度分秒连写的(double型)角度化为弧度值
double dms_rad(double a)
{
	//提取角度值的符号
	double sign=(a<0.0) ? -1.0 : 1.0;
	a=fabs(a);
	
	//提取角度值的整度
	int d=(int)((a+0.00001)/10000.0);
	a=a-d*10000.0;
	if(a<0.0){ d=d-1; a=a+10000;}
	
	//提取角度值的整分及秒值
	int m=(int)((a+0.00001)/100.0);
	a=a-m*100;
	if(a<0.0){ m=m-1; a=a+100.0;}
	
	a=sign*(d*3600.0+m*60.0+a)/206264.806247096363;
	
	return a;
}


//////////////////////////////////////////////////////////////////////////
//  将角度的弧度值化为度分秒连写的角度（double 型） 
double rad_dms(double a)
{
	a=a*206264.806247096363;
	
	double sign=(a<0.0) ? -1.0 : 1.0;
	a=fabs(a);
	
	int d=(int)(a/3600.0+0.0000001);
	a=a-d*3600.0;
	
	if(a<0.0){ d=d-1; a=a+3600.0; }
	
	int m=(int)(a/60.0+0.0001);
	a=a-m*60.0;
	if(a<0.0){ m=m-1; a=a+60.0; }
	
	a=d*10000.0+m*100.0+a;
	
	return a*sign;
}



//////////////////////////////////////////////////////////////////////////
//    权逆阵传播计算
void Calculate_BQBT(double B[],double Q[],int r,int n,double N[])
{
	for(int i=0;i<r;i++)
		for(int j=0;j<=i;j++)
		{
			double nij=0.0;
			for(int k=0;k<n;k++)
				for(int s=0;s<n;s++)
					nij+=B[i*n+k]*Q[ij(k,s)]*B[j*n+s];
				N[ij(i,j)]+=nij;
		}
}




//////////////////////////////////////////////////////////////////////////
//    向文件输出数组
void PrintM(FILE *fp,double A[],int size, int t,char* fmt,
			char* title,bool IsLabel)
{
	if(title)fprintf(fp,"\n %s: ",title);
	int j=0;
	for(int i=0;i<size;i++)
	{
		if(i%t==0)
		{
			j++;
			if(IsLabel)fprintf(fp,"\n%3d ",j);
			else fprintf(fp,"\n    ");
		}
		fprintf(fp,fmt,A[i]);
	}	
	fprintf(fp,"\n");
}



//////////////////////////////////////////////////////////////////////////
//  向文件输出对称矩阵（下三角元素）
void  PrintM2(FILE* fp, double M[], int n, int t,const char *fmt,
			  const char* title,bool IsLabel)
{
	if(title)fprintf(fp,"\n %s: ",title);
	
	int index=0;
	for(int i=0;i<n;i++)
	{
		if(IsLabel)fprintf(fp,"\n%3d ",i+1);
		else fprintf(fp,"\n    ");
		for(int j=0;j<=i;j++)
		{
			if(j>0 && j%t==0)fprintf(fp,"\n    ");
			fprintf(fp,fmt,M[index++]);
		}
	}
	fprintf(fp,"\n");	
}


//////////////////////////////////////////////////////////////////////////
//  向文件输出线性方程组
void  PrintEquation(FILE* fp, double A[], double b[],
					int n, int t, char *fmt, char* title)
{
	if(title)fprintf(fp,"\n %s: ",title);
	for(int i=0;i<n;i++)
	{
		fprintf(fp,"\n%3d ",i+1);
		for(int j=0;j<t;j++)
		{
			fprintf(fp,fmt,A[i*t+j]);
		}
		fprintf(fp,fmt,b[i]);
	}
}


//////////////////////////////////////////////////////////////////////////
//  对称正定矩阵求逆(仅存下三角元素)
bool inverse(double a[],int n)
{
    double *a0=new double[n];
    for(int k=0;k<n;k++)
	{ 
		double a00=a[0];
		if(a00+1.0==1.0)
		{
			delete []a0; 
			return false;
		}
		for(int i=1;i<n;i++)
		{
			double ai0 = a[i*(i+1)/2];
			
			if(i<=n-k-1)a0[i]= -ai0/a00;
			else        a0[i]=  ai0/a00;

			for(int j=1;j<=i;j++)
			{
				a[(i-1)*i/2+j-1]=a[i*(i+1)/2+j]+ai0*a0[j];
			} 
		}
		for(int i=1;i<n;i++)
		{
			a[(n-1)*n/2+i-1]=a0[i];
		}
		a[n*(n+1)/2-1]=1.0/a00;
	} 
	delete []a0;
	return true;
}


#define IGG1  0
#define IGG3  1
#define Huber 2
double Wi(int fname, double v, double k0, double k1)
{
	double a;
	switch(fname)
	{
	case IGG1://IGG1函数
		v=fabs(v);
		if(v<=k0)return 1.0;
		if(v>k1)return 0.0;
		return k0/v;
	case IGG3://IGG3函数
		v=fabs(v);
		if(v<=k0)return 1.0;
		if(v>k1)return 0.0;
		a=(k1-v)/(k1-k0);
		return k0/v*a*a;
	case Huber: // Huber函数
		v=fabs(v);
		if(v<=k0)return 1.0;
		return k0/v;
	default:
		MyBreak("等价权函数名称错误！");
	return 1.0;
	}
}



//////////////////////////////////////////////////////////////////////////
// 中位数计算：将数组由大至小排序，返回值为中位数
double Median(double pp[],int n,bool IsAbs)
{
	double *p=new double [n];
	if(IsAbs)
	{
		for(int i=0;i<n;i++)p[i]=fabs(pp[i]);
	}
	else
	{
		for(int i=0;i<n;i++)p[i]=pp[i];
	}
	
	int k=n/2;
	while(k>0)
	{  
		for(int j=k;j<=n-1;j++)
		{  
			double t=p[j]; 
			int i=j-k;
			
			while( (i>=0) && (p[i]>t) )
			{
				p[i+k]=p[i];	
				i=i-k;  
			}
			p[i+k]=t;
		}
		k=k/2;
	}
	
	double mean = (n%2==1)? p[n/2] : (p[n/2]+p[n/2-1])/2.0;
	delete []p;
	
	return  mean;
	
}

struct FourPara
{
	double dx;
	double dy;
	double s;//比例尺
	double r;//旋转角
};
struct XYPoint
{
	std::string name;
	double x;
	double y;
};
//左手系，逆时针旋转
FourPara LS_PLCC(const std::vector<XYPoint>& source, const std::vector<XYPoint>& target, const int Num, const bool b_scale = false)
{
	FourPara fourpara;
	double X0(0), Y0(0), S0(1), R0(0);//四参数初始值
	Eigen::MatrixXd A(2 * Num, 4);
	if (b_scale)
	{
		A.resize(2 * Num, 3);
	}
	A.setZero();//初始化矩阵全部为0
	Eigen::MatrixXd P(2 * Num, 2 * Num);
	P.setIdentity();//设置为单位矩阵
	Eigen::VectorXd L(2 * Num);

	//计算四参数初值
	double max = 0;//最远距离
	int n_mx = 0;//最远距离索引
	for (size_t i = 1; i < source.size(); i++)
	{
		double temlen = sqrt(pow((source[i].y - source[0].y), 2) + pow((source[i].x - source[0].x), 2));
		if (temlen > max)
		{
			max = temlen;
			n_mx = i;
		}
	}
	double dX1 = source[n_mx].x - source[0].x;
	double dY1 = source[n_mx].y - source[0].y;
	double dX2 = target[n_mx].x - target[0].x;
	double dY2 = target[n_mx].y - target[0].y;

	double dR1 = sqrt(dX1 * dX1 + dY1 * dY1);
	dX1 = dX1 / dR1;
	dY1 = dY1 / dR1;
	double dR2 = sqrt(dX2 * dX2 + dY2 * dY2);
	dX2 = dX2 / dR2;
	dY2 = dY2 / dR2;

	double RR = dX1 * dY2 - dX2 * dY1;
	if (abs(RR) < 1e-10)
	{
		if (abs(dX1 - dX2) < 1e-8 && abs(dY1 - dY2) < 1e-8)
		{
			R0 = 0;
		}
		else
		{
			R0 = PI;
		}
	}
	if (RR < 0)
	{
		R0 = asin(RR) + 2 * PI;
	}
	else
	{
		R0 = asin(RR);
	}

	X0 = Y0 = 0;
	S0 = 1;
	Eigen::Vector4d ds(1, 1, 1, 1);
	int tatolnum = 0;//循环次数
	while (!(fabs(ds(0)) < 0.000001 && fabs(ds(1)) < 0.000001 && fabs(ds(2)) < 0.000001 && fabs(ds(3)) < 0.00000001))
	{
		//设置系数矩阵
		for (size_t i = 0; i < Num; i++)
		{
			A(2 * i, 0) = 1;
			A(2 * i, 1) = 0;
			A(2 * i, 2) = S0 * (-sin(R0) * source[i].x - cos(R0) * source[i].y);
			if (!b_scale)//固定比例尺
			{
				A(2 * i, 3) = cos(R0) * source[i].x - sin(R0) * source[i].y;
			}
			A(2 * i + 1, 0) = 0;
			A(2 * i + 1, 1) = 1;
			A(2 * i + 1, 2) = S0 * (cos(R0) * source[i].x - sin(R0) * source[i].y);
			if (!b_scale)
			{
				A(2 * i + 1, 3) = sin(R0) * source[i].x + cos(R0) * source[i].y;
			}

			L(2 * i) = target[i].x - (X0 + S0 * (cos(R0) * source[i].x - sin(R0) * source[i].y));
			L(2 * i + 1) = target[i].y - (Y0 + S0 * (sin(R0) * source[i].x + cos(R0) * source[i].y));
		}
		//计算atpa
		Eigen::MatrixXd ATPA = A.transpose() * P * A;
		if (abs(ATPA.determinant()) < 1e-10)
		{
			std::cout << "ATPA不可逆" << std::endl;
		}
		Eigen::VectorXd ds1 = ATPA.inverse() * (A.transpose() * P * L);
		for (size_t i = 0; i < ds1.size(); i++)
		{
			ds(i) = ds1(i);
		}
		if (b_scale)ds(3) = 0;
		X0 += ds(0);
		Y0 += ds(1);
		R0 += ds(2);
		if (!b_scale)
		{
			S0 += ds(3);
		}

		if (R0 > 2 * PI)
		{
			R0 -= 2 * PI;
		}
		if (R0 < 0)
		{
			R0 += 2 * PI;
		}

		fourpara.dx = X0;
		fourpara.dy = Y0;
		fourpara.s = S0;
		fourpara.r = R0;
		tatolnum++;
		if (tatolnum > 100)
		{
			std::cout << "超过100次" << std::endl;
			break;
		}
	}
	return fourpara;
}

XYPoint Trans_PLPoint(const XYPoint& p, const FourPara fp)
{
	XYPoint p1;
	p1.x = fp.dx + fp.s * (cos(fp.r) * p.x - sin(fp.r) * p.y);
	p1.y = fp.dy + fp.s * (sin(fp.r) * p.x + cos(fp.r) * p.y);
	return p1;
}
// 计算四参数转换 (左手系)
void calculateTransformationParameters(const int n, const double* srcPoints, const double* dstPoints, double& dx, double& dy, double& scale, double& rotation) {

	XYPoint p1;
	std::vector<XYPoint> sp;
	std::vector<XYPoint> tp;
	for (size_t i = 0; i < n; i++)
	{
		p1.x = srcPoints[2 * i];
		p1.y = srcPoints[2 * i + 1];
		sp.push_back(p1);
		p1.x = dstPoints[2 * i];
		p1.y = dstPoints[2 * i + 1];
		tp.push_back(p1);
	}
	FourPara fp = LS_PLCC(sp, tp, n, true);
	scale = fp.s;
	rotation = fp.r;
	dx = fp.dx;
	dy = fp.dy;
}

// 使用四参数进行坐标转换 
void transformPoint(const double x, const double y, double& rex, double& rey, double dx, double dy, double scale, double rotation) {

	rex = dx + scale * (x * cos(rotation) - y * sin(rotation));
	rey = dy + scale * (x * sin(rotation) + y * cos(rotation));
}