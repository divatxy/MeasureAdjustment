#include "../MeasureAdjustment/MeasureAdjustment.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include <math.h>
#include <iomanip>
//#include <numbers>
using namespace std;
#define PI 3.14159265358979323846
//��ȡ�ļ�������۲����ݣ���SVD����Ϊˮƽ�ǡ�ƽ��͸߲�
struct OrientData
{
    std::string staName;
    std::vector<string> tarName;
    vector<double> oriendata;
    vector<double> v_orient;//������
    vector<double> mu_orient;//�����
};
struct LengthData
{
    string staName;
    string tarName;
    double lenData;
    double v_len;//������
    double mu_len;//�����
};
struct PXYZ
{
    string name;
    double X;
    double Y;
    double Z;
    double mx;
    double my;
    double mz;
    double mM;
};

double rad_deg(const double rad)
{
    return rad * 360 / PI;
}
double deg_rad(const double deg)
{
    return deg * PI / 360.0;
}
//////////////////////////////////////////////////////////////////////////
//  ���ȷ�����д��(double��)�ǶȻ�Ϊ����ֵ
double dms_rad(double a)
{
    //��ȡ�Ƕ�ֵ�ķ���
    double sign = (a < 0.0) ? -1.0 : 1.0;
    a = fabs(a);

    //��ȡ�Ƕ�ֵ������
    int d = (int)((a + 0.00001) / 10000.0);
    a = a - d * 10000.0;
    if (a < 0.0) { d = d - 1; a = a + 10000; }

    //��ȡ�Ƕ�ֵ�����ּ���ֵ
    int m = (int)((a + 0.00001) / 100.0);
    a = a - m * 100;
    if (a < 0.0) { m = m - 1; a = a + 100.0; }

    a = sign * (d * 3600.0 + m * 60.0 + a) / 206264.806247096363;

    return a;
}


//////////////////////////////////////////////////////////////////////////
//  ���ǶȵĻ���ֵ��Ϊ�ȷ�����д�ĽǶȣ�double �ͣ� 
double rad_dms(double a)
{
    a = a * 206264.806247096363;

    double sign = (a < 0.0) ? -1.0 : 1.0;
    a = fabs(a);

    int d = (int)(a / 3600.0 + 0.0000001);
    a = a - d * 3600.0;

    if (a < 0.0) { d = d - 1; a = a + 3600.0; }

    int m = (int)(a / 60.0 + 0.0001);
    a = a - m * 60.0;
    if (a < 0.0) { m = m - 1; a = a + 60.0; }

    a = d * 10000.0 + m * 100.0 + a;

    return a * sign;
}


void planeAD(const int knNumber,vector<PXYZ> &knPoint,vector<OrientData> &orienVec,vector<LengthData> &lenVec )
{
    int Pnumber(0), knPnumber(0), Lnumber(0), Nnumber(0), Snumber(0), Tnumber(0);
    double ma(1), ms1(0.005), ms2(1.0e-6), mt(0);
    Pnumber = knPoint.size();
    knPnumber = knNumber;
    Lnumber = orienVec.size();
    for (auto& item : orienVec)
    {
        Nnumber += item.oriendata.size();
    }
    Snumber = lenVec.size();

    char** pname = new char* [Pnumber];//����
    for (size_t i = 0; i < Pnumber; i++)
    {        
        pname[i] = (char*)knPoint[i].name.c_str();
    }
    //�����꼰���Ӧ����
    double* XY = new double[2 * Pnumber];
    double* mx = new double[Pnumber];
    double* my = new double[Pnumber];
    double* M = new double[Pnumber];

    int* dir0{}; int* dir1{}; int* dir2{}; int* S_dir1{}; int* S_dir2{}; int* T_dir1{}; int* T_dir2{};
    double* L{}; double* V{}; double* R_T{}; double* R_mt{}; double* R_S{}; double* R_ms{};
    double* S_L{}; double* S_V{}; double* S_T{}; double* S_mt{}; double* S_S{}; double* S_ms{};
    double* T_L{}; double* T_V{}; double* T_T{}; double* T_mt{}; double* T_S{}; double* T_ms{};

    //����۲�ֵ����
    if (Lnumber > 0)
    {
        dir0 = new int[Lnumber + 1]; //�����׷���۲�ֵ�����
        dir1 = new int[Lnumber];   //��վ���
        dir2 = new int[Nnumber];   //��׼���
        L = new double[Nnumber];   //����ֵ
        V = new double[Nnumber];   //�ֲ�
        R_T = new double[Nnumber];
        R_mt = new double[Nnumber];
        R_S = new double[Nnumber];
        R_ms = new double[Nnumber];
    }
    if (Snumber > 0)//Ϊ�߳��۲�ֵ���������ڴ�
    {
        S_dir1 = new int[Snumber];  //��վ���
        S_dir2 = new int[Snumber];  //��׼���
        S_L = new double[Snumber];  //�߳��۲�ֵ
        S_V = new double[Snumber];  //�߳��в�
        S_T = new double[Snumber];
        S_mt = new double[Snumber];
        S_S = new double[Snumber];
        S_ms = new double[Snumber];
    }
    if (Tnumber > 0)//Ϊ��λ�ǹ۲�ֵ���������ڴ�
    {
        T_dir1 = new int[Tnumber];  //��վ���
        T_dir2 = new int[Tnumber];  //��׼���
        T_L = new double[Tnumber];  //��λ�ǹ۲�ֵ
        T_V = new double[Tnumber];  //��λ�ǲв�
        T_T = new double[Tnumber];
        T_mt = new double[Tnumber];
        T_S = new double[Tnumber];
        T_ms = new double[Tnumber];
    }

    for (size_t i = 0; i < knNumber; i++)//��֪��
    {
        XY[2 * i] = knPoint[i].X;
        XY[2 * i + 1] = knPoint[i].Y;
    }
    
    //��ȡ����ֵ
    dir0[0] = 0;
    for (size_t i = 0; i < orienVec.size(); i++)
    {
        for (size_t i1 = 0; i1 < knPoint.size(); i1++)
        {
            if (knPoint[i1].name==orienVec[i].staName)
            {
                dir1[i] = i1; 
                break;
            }
        }
        dir0[i + 1] = dir0[i] + orienVec[i].tarName.size();
        for (size_t i2 = 0; i2 < orienVec[i].tarName.size(); i2++)
        {
            for (size_t i3 = 0; i3 < knPoint.size(); i3++)
            {
                if (knPoint[i3].name == orienVec[i].tarName[i2])
                {
                    dir2[dir0[i]+i2] = i3;//��׼���
                    break;
                }
            }
            L[dir0[i] + i2] = orienVec[i].oriendata[i2] * PI / 360.0;//�Ƕ�תΪ����
        }
    }
        
    //��ȡ�߳�
    for (size_t i = 0; i < lenVec.size(); i++)
    {
        for (size_t i1 = 0; i1 < knPoint.size(); i1++)
        {
            if (knPoint[i1].name == lenVec[i].staName)
            {
                S_dir1[i] = i1;//��˵���
            }
            if (knPoint[i1].name == lenVec[i].tarName)
            {
                S_dir2[i] = i1;//�Ҷ˵���
            }
        }
        S_L[i] = lenVec[i].lenData;
    }

    double mu = 0;//ƽ���ܾ���
    PlaneNetAD1(Pnumber, knPnumber, Lnumber, Nnumber, Snumber, Tnumber, ma, ms1, ms2, mt, pname,
        dir0, dir1, dir2, L, V, R_T, R_mt, R_S, R_ms,
        S_dir1, S_dir2, S_L, S_V, S_T, S_mt, S_S, S_ms,
        T_dir1, T_dir2, T_L, T_V, T_T, T_mt, T_S, T_ms,
        XY, mx, my, M, mu);

    std::cout << "��С����ƽ��:��=" << mu << std::endl;
    std::cout << "------����ƽ��侫��-----" << std::endl;
    std::cout << "No.  P   X   Y   mx   my  M" << std::endl;
    for (size_t i = 0; i < Pnumber; i++)
    {
        std::cout << i << "\t" << pname[i] << "\t" << XY[2 * i] << "\t" << XY[2 * i + 1] << "\t" << mx[i] << "\t" << my[i] << "\t" << M[i] << std::endl;
    }
    std::cout << "\n" << std::endl;
    std::cout << "------����۲�ֵƽ��ɹ�-----" << std::endl;
    std::cout << "P1 \t P2 \t L \t V \t T \t mT \t S \t ms" << std::endl;
    for (size_t i = 0; i < Lnumber; i++)
    {
        //std::cout << "\n" << std::endl;
        //std::cout << pname[dir1[i]];
        for (size_t j = dir0[i]; j < dir0[i + 1]; j++)
        {
            if (j == dir0[i])
            {
                std::cout << "\n" << pname[dir1[i]];

            }
            else
            {
                //std::cout << "\n";
            }
            std::cout << "\t" << pname[dir2[j]] << "\t" << rad_deg(L[j]) << "\t" << V[j] << "\t" << rad_deg(R_T[j]) << "\t" << R_mt[j] << "\t" << R_S[j] << "\t" << R_ms[j] << std::endl;
        }
    }

    std::cout << "\n" << std::endl;
    std::cout << "------�߳��۲�ֵƽ��ɹ�-----" << std::endl;
    std::cout << "P1 \t P2 \t S \t Vs \t T \t mT \t S+Vs \t ms" << std::endl;
    for (size_t i = 0; i < Snumber; i++)
    {
        std::cout << pname[S_dir1[i]] << "\t" << pname[S_dir2[i]] << "\t" << S_L[i] << "\t" << S_V[i] << "\t" << rad_deg(S_T[i]) << "\t" << S_mt[i] << "\t" << S_S[i] << "\t" << S_ms[i] << std::endl;
    }

    

}



void threeAD()
{
    std::ifstream file("..\\ThreeNetAD\\����վ��������20240811-1.txt");
    // ����ļ��Ƿ�ɹ���
    if (!file.is_open()) {
        std::cerr << "�޷����ļ� " << std::endl;
        return; // ����ļ���ʧ�ܣ��˳�����
    }
    vector<string> allPointName;
    string sta = "";
    string tempstastr("");
    string temptarstr("");
    double Hv(0), V(0), S(0);
    double H(0), S_H(0);
    double _hv(0);//ˮƽ�͹���Ĳ�ֵ
    std::string line;
    std::getline(file, line);//������һ��
    OrientData od;
    vector<OrientData> vod;//����洢����۲�ֵ
    vector<LengthData> vld;//�洢�߳��۲�ֵ
    while (std::getline(file, line))
	{
        LengthData ld;//�߳�����
		if (line == "")
		{
			break;
		}
		std::istringstream iss(line);
		
		iss >> tempstastr >> temptarstr >> Hv >> V >> S;
		bool is_staIn(true), is_tarIn(true);

		//----------------�����洢-----------------
		for (auto& i : allPointName)
		{
			if (i == tempstastr && is_staIn)
			{
				is_staIn = false;
			}
			if (i == temptarstr && is_tarIn)
			{
				is_tarIn = false;
			}
			if (!is_staIn && !is_tarIn)
			{
				break;
			}
		}
		if (is_staIn)
		{
			allPointName.push_back(tempstastr);
		}
		if (is_tarIn)
		{
			allPointName.push_back(temptarstr);
		}
		//----------------�����洢����-----------------

		if (sta != tempstastr)//һ��Ŀ�ʼ
		{
			if (od.staName != "")
			{
				vod.push_back(od);
			}
			od.staName.clear();
			od.tarName.clear();
			od.oriendata.clear();
			sta = tempstastr;
			_hv = Hv;
			Hv = 0;
			od.staName = sta;
			od.tarName.push_back(temptarstr);
			od.oriendata.push_back(Hv);
            ld.staName = tempstastr;
            ld.tarName = temptarstr;
            ld.lenData = sin(V*PI/360) * S/1000.0;//תΪ��
            vld.push_back(ld);
		}
		else
		{
			Hv = Hv - _hv;
			if (Hv < 0)
			{
				Hv = 360 + Hv;
			}
			else if (Hv > 360)
			{
				Hv = Hv - 360;
			}
			od.tarName.push_back(temptarstr);
			od.oriendata.push_back(Hv);
            ld.staName = tempstastr;
            ld.tarName = temptarstr;
            ld.lenData = sin(V * PI / 360) * S/1000.0;
            vld.push_back(ld);
		}
	}
	//�����һ��۲�ֵ������������
	vod.push_back(od);

	//-----------��ȡ��ֵ֪-------------
	std::ifstream file1("..\\ThreeNetAD\\��֪�����꣨�߳�Ϊˮ׼�̣߳�20240811(1).txt");
	// ����ļ��Ƿ�ɹ���
	if (!file1.is_open()) {
		std::cerr << "�޷����ļ� " << std::endl;
		return; // ����ļ���ʧ�ܣ��˳�����
	}
	std::getline(file1, line);//������һ��
	vector<PXYZ> P_sta;
	double X(0), Y(0);
	/*while (std::getline(file1, line))
	{
		if (line == "")
		{
			break;
		}
		PXYZ xyz;
		std::istringstream iss(line);
		iss >> sta >> X >> Y >> H;
		xyz.name = sta;
		xyz.X = X;
		xyz.Y = Y;
		xyz.Z = H;
		P_sta.push_back(xyz);
		allPointName.erase(
			std::remove_if(allPointName.begin(), allPointName.end(), [sta](const std::string& s) {return s == sta; }),
			allPointName.end()
		);
	}*/
	int knNum = P_sta.size();
	for (auto& item : allPointName)
	{
		PXYZ p;
		p.name = item;
		P_sta.push_back(p);
	}
    
    planeAD(knNum, P_sta, vod, vld);

	//��ʾvod����---
	/*for (auto&e : vod)
	{
		cout << e.staName << endl;
		for (size_t i = 0; i < e.tarName.size(); i++)
		{
			cout << "\t" << e.tarName[i] << "\t" << e.oriendata[i] << endl;
		}
	}*/
	//----��ʾ����---
    file.close();
    file1.close();



    return;

    //д�����ݵ��ı��ĵ�
    std::ofstream file2("..\\TempOutData1.txt");
    // ����ļ��Ƿ�ɹ���
    if (!file2.is_open()) {
        std::cerr << "�޷����ļ� " << std::endl;
        return; // ����ļ���ʧ�ܣ��˳�����
    }
    int nNum = 0;
    for (auto& it : vod)
    {
        for (auto& ite : it.oriendata)
        {
            nNum++;
        }
    }
	file2 << P_sta.size() << "," << knNum << "," << vod.size() << "," << nNum << "," << vld.size() << "," << 0 << endl;
    file2 << 1.0 << "," << 0.005 << "," << 1.0e-6 << "," << 0 <<  endl;
	file2 << std::fixed << std::endl;
	for (size_t i = 0; i < knNum; i++)
	{
		file2 << P_sta[i].name << ",";
		file2 << std::setprecision(4) << P_sta[i].X << "," << P_sta[i].Y << std::endl;
	}
	//file2   <<   std::endl;
	for (auto& it : vod)
	{
		file2 << std::endl;//����
		file2 << it.staName << std::endl;
		for (size_t i = 0; i < it.tarName.size(); i++)
		{
			file2 << it.tarName[i] << ",L," << std::setprecision(2) << rad_dms(deg_rad(it.oriendata[i])) << endl;
		}
        for (auto& ite : vld)
        {
            if (ite.staName==it.staName)
            {
                file2<<ite.tarName<<",S,"<< std::setprecision(4) << ite.lenData << endl;
            }
        }
	}
	//file2 << std::endl;//����
	//for (auto& it : vld)
	//{
	//	file2 << it.staName << "   " << it.tarName << "   " << std::setprecision(4) << it.lenData << endl;
	//}

	file2.close();
}


int main()
{
	threeAD();
    std::cout << "Hello World!\n";
}

// ���г���: Ctrl + F5 ����� >����ʼִ��(������)���˵�
// ���Գ���: F5 ����� >����ʼ���ԡ��˵�

// ����ʹ�ü���: 
//   1. ʹ�ý��������Դ�������������/�����ļ�
//   2. ʹ���Ŷ���Դ�������������ӵ�Դ�������
//   3. ʹ��������ڲ鿴���������������Ϣ
//   4. ʹ�ô����б��ڲ鿴����
//   5. ת������Ŀ��>���������Դ����µĴ����ļ�����ת������Ŀ��>�����������Խ����д����ļ���ӵ���Ŀ
//   6. ��������Ҫ�ٴδ򿪴���Ŀ����ת�����ļ���>���򿪡�>����Ŀ����ѡ�� .sln �ļ�
