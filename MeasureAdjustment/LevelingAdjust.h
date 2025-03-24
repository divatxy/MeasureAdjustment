// LevelingAdjust.h: interface for the CLevelingAdjust class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LEVELINGADJUST_H__553D83A6_1E43_4AB0_9C3E_8007847A4AA0__INCLUDED_)
#define AFX_LEVELINGADJUST_H__553D83A6_1E43_4AB0_9C3E_8007847A4AA0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "pch.h"

#include "stdio.h"
#include <string>
#include <vector>

class CLevelingAdjust  
{
public:
	CLevelingAdjust();
	virtual ~CLevelingAdjust();

	int m_Lnumber;      //�߲�����
	int m_Pnumber;      //�ܵ���
	int m_knPnumber;    //��֪����
	int m_StablePnumber; //���ȵ���
	double m_pvv;        //[pvv]
	FILE   *resultfp;    //�ļ�ָ��,���������

	double  m_Sigma;  //��ǰ��λȨ�����ֲ�̽�⡢ �պϲ�����ã��ڴֲ�̽�⡢�պϲ���һ��Ϊ0.001����������������С����һ��Ϊ0.005
	double  Alhpa;     // ����ˮƽ���ֲ�̽������ã�
	
	int *StartP;     //�߲�����
	int *EndP;       //�߲��յ��
	char   **Pname;   //������ַ����
	double  *L;       //�۲�ֵ����
	double  *Height;  //�߳�ֵ����
	double  *P;       //�۲�ֵ��Ȩ

	double  *ATPA,*ATPL; //������ϵ��������������
	double  *dX;      //����ƽ��ֵ���̸߳�������   
	double  *V;       //�в�����
	double  m_mu;      //���λȨ�����
	
	int *IsStable; //�Ƿ�Ϊ���ȵ��
		
	void   Inputdata(char *datafile);//����ԭʼ���ݺ���
	void   Printdata(); //���ԭʼ���ݺ���
	int    GetStationNumber(char *name); //��ȡ��ź���
	void   ca_H0();       //���Ƹ̼߳��㺯��
	void   ca_ATPA();  //��������ɺ���
	void   ca_dX();    //ƽ��ֵ���㺯��
	void   PrintResult();  //���ȹ�����ƽ��ֵ�������
	double ca_V();     //�в���㺯��	

    void   LS_Adjustment();//��С����ƽ���
	void   Quasi_Stable( char *file);//����ƽ���
    void   FreeNetAdjust();//������ƽ���
		
	void  FindShortPath(int p1,int exclude, int root[],
						  double diff[],double S[]);    //���·������
	int  LineClosure(std::vector<std::string>& line_name, std::vector<double>& line_L, std::vector<double>& line_w, std::vector<double>& line_limit); // ·�߱պϲ����
	int  LoopClosure(std::vector<std::string>& loop_name, std::vector<double>& loop_L, std::vector<double>& loop_w, std::vector<double>& loop_limit); //���պϲ����
	

	void DataSnooping(double arfa);//�ֲ�̽��
		
};

#endif // !defined(AFX_LEVELINGADJUST_H__553D83A6_1E43_4AB0_9C3E_8007847A4AA0__INCLUDED_)





















