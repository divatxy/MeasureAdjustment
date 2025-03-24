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
//读取文件并整理观测数据，将SVD整理为水平角、平距和高差
struct OrientData
{
    std::string staName;
    std::vector<string> tarName;
    vector<double> oriendata;
    vector<double> v_orient;//改正数
    vector<double> mu_orient;//中误差
};
struct LengthData
{
    string staName;
    string tarName;
    double lenData;
    double v_len;//改正数
    double mu_len;//中误差
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
//  将度分秒连写的(double型)角度化为弧度值
double dms_rad(double a)
{
    //提取角度值的符号
    double sign = (a < 0.0) ? -1.0 : 1.0;
    a = fabs(a);

    //提取角度值的整度
    int d = (int)((a + 0.00001) / 10000.0);
    a = a - d * 10000.0;
    if (a < 0.0) { d = d - 1; a = a + 10000; }

    //提取角度值的整分及秒值
    int m = (int)((a + 0.00001) / 100.0);
    a = a - m * 100;
    if (a < 0.0) { m = m - 1; a = a + 100.0; }

    a = sign * (d * 3600.0 + m * 60.0 + a) / 206264.806247096363;

    return a;
}


//////////////////////////////////////////////////////////////////////////
//  将角度的弧度值化为度分秒连写的角度（double 型） 
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

    char** pname = new char* [Pnumber];//点名
    for (size_t i = 0; i < Pnumber; i++)
    {        
        pname[i] = (char*)knPoint[i].name.c_str();
    }
    //点坐标及其对应精度
    double* XY = new double[2 * Pnumber];
    double* mx = new double[Pnumber];
    double* my = new double[Pnumber];
    double* M = new double[Pnumber];

    int* dir0{}; int* dir1{}; int* dir2{}; int* S_dir1{}; int* S_dir2{}; int* T_dir1{}; int* T_dir2{};
    double* L{}; double* V{}; double* R_T{}; double* R_mt{}; double* R_S{}; double* R_ms{};
    double* S_L{}; double* S_V{}; double* S_T{}; double* S_mt{}; double* S_S{}; double* S_ms{};
    double* T_L{}; double* T_V{}; double* T_T{}; double* T_mt{}; double* T_S{}; double* T_ms{};

    //方向观测值数组
    if (Lnumber > 0)
    {
        dir0 = new int[Lnumber + 1]; //各组首方向观测值的序号
        dir1 = new int[Lnumber];   //测站点号
        dir2 = new int[Nnumber];   //照准点号
        L = new double[Nnumber];   //方向值
        V = new double[Nnumber];   //粗差
        R_T = new double[Nnumber];
        R_mt = new double[Nnumber];
        R_S = new double[Nnumber];
        R_ms = new double[Nnumber];
    }
    if (Snumber > 0)//为边长观测值数组申请内存
    {
        S_dir1 = new int[Snumber];  //测站点号
        S_dir2 = new int[Snumber];  //照准点号
        S_L = new double[Snumber];  //边长观测值
        S_V = new double[Snumber];  //边长残差
        S_T = new double[Snumber];
        S_mt = new double[Snumber];
        S_S = new double[Snumber];
        S_ms = new double[Snumber];
    }
    if (Tnumber > 0)//为方位角观测值数组申请内存
    {
        T_dir1 = new int[Tnumber];  //测站点号
        T_dir2 = new int[Tnumber];  //照准点号
        T_L = new double[Tnumber];  //方位角观测值
        T_V = new double[Tnumber];  //方位角残差
        T_T = new double[Tnumber];
        T_mt = new double[Tnumber];
        T_S = new double[Tnumber];
        T_ms = new double[Tnumber];
    }

    for (size_t i = 0; i < knNumber; i++)//已知点
    {
        XY[2 * i] = knPoint[i].X;
        XY[2 * i + 1] = knPoint[i].Y;
    }
    
    //读取方向值
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
                    dir2[dir0[i]+i2] = i3;//照准点号
                    break;
                }
            }
            L[dir0[i] + i2] = orienVec[i].oriendata[i2] * PI / 360.0;//角度转为弧度
        }
    }
        
    //读取边长
    for (size_t i = 0; i < lenVec.size(); i++)
    {
        for (size_t i1 = 0; i1 < knPoint.size(); i1++)
        {
            if (knPoint[i1].name == lenVec[i].staName)
            {
                S_dir1[i] = i1;//左端点名
            }
            if (knPoint[i1].name == lenVec[i].tarName)
            {
                S_dir2[i] = i1;//右端点名
            }
        }
        S_L[i] = lenVec[i].lenData;
    }

    double mu = 0;//平差总精度
    PlaneNetAD1(Pnumber, knPnumber, Lnumber, Nnumber, Snumber, Tnumber, ma, ms1, ms2, mt, pname,
        dir0, dir1, dir2, L, V, R_T, R_mt, R_S, R_ms,
        S_dir1, S_dir2, S_L, S_V, S_T, S_mt, S_S, S_ms,
        T_dir1, T_dir2, T_L, T_V, T_T, T_mt, T_S, T_ms,
        XY, mx, my, M, mu);

    std::cout << "最小二乘平差:μ=" << mu << std::endl;
    std::cout << "------坐标平差及其精度-----" << std::endl;
    std::cout << "No.  P   X   Y   mx   my  M" << std::endl;
    for (size_t i = 0; i < Pnumber; i++)
    {
        std::cout << i << "\t" << pname[i] << "\t" << XY[2 * i] << "\t" << XY[2 * i + 1] << "\t" << mx[i] << "\t" << my[i] << "\t" << M[i] << std::endl;
    }
    std::cout << "\n" << std::endl;
    std::cout << "------方向观测值平差成果-----" << std::endl;
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
    std::cout << "------边长观测值平差成果-----" << std::endl;
    std::cout << "P1 \t P2 \t S \t Vs \t T \t mT \t S+Vs \t ms" << std::endl;
    for (size_t i = 0; i < Snumber; i++)
    {
        std::cout << pname[S_dir1[i]] << "\t" << pname[S_dir2[i]] << "\t" << S_L[i] << "\t" << S_V[i] << "\t" << rad_deg(S_T[i]) << "\t" << S_mt[i] << "\t" << S_S[i] << "\t" << S_ms[i] << std::endl;
    }

    

}



void threeAD()
{
    std::ifstream file("..\\ThreeNetAD\\各测站测量数据20240811-1.txt");
    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << std::endl;
        return; // 如果文件打开失败，退出程序
    }
    vector<string> allPointName;
    string sta = "";
    string tempstastr("");
    string temptarstr("");
    double Hv(0), V(0), S(0);
    double H(0), S_H(0);
    double _hv(0);//水平就归零的差值
    std::string line;
    std::getline(file, line);//跳过第一行
    OrientData od;
    vector<OrientData> vod;//安组存储方向观测值
    vector<LengthData> vld;//存储边长观测值
    while (std::getline(file, line))
	{
        LengthData ld;//边长数据
		if (line == "")
		{
			break;
		}
		std::istringstream iss(line);
		
		iss >> tempstastr >> temptarstr >> Hv >> V >> S;
		bool is_staIn(true), is_tarIn(true);

		//----------------点名存储-----------------
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
		//----------------点名存储结束-----------------

		if (sta != tempstastr)//一组的开始
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
            ld.lenData = sin(V*PI/360) * S/1000.0;//转为米
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
	//将最后一组观测值加入分组变量中
	vod.push_back(od);

	//-----------读取已知值-------------
	std::ifstream file1("..\\ThreeNetAD\\已知点坐标（高程为水准高程）20240811(1).txt");
	// 检查文件是否成功打开
	if (!file1.is_open()) {
		std::cerr << "无法打开文件 " << std::endl;
		return; // 如果文件打开失败，退出程序
	}
	std::getline(file1, line);//跳过第一行
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

	//显示vod数据---
	/*for (auto&e : vod)
	{
		cout << e.staName << endl;
		for (size_t i = 0; i < e.tarName.size(); i++)
		{
			cout << "\t" << e.tarName[i] << "\t" << e.oriendata[i] << endl;
		}
	}*/
	//----显示结束---
    file.close();
    file1.close();



    return;

    //写入数据到文本文档
    std::ofstream file2("..\\TempOutData1.txt");
    // 检查文件是否成功打开
    if (!file2.is_open()) {
        std::cerr << "无法打开文件 " << std::endl;
        return; // 如果文件打开失败，退出程序
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
		file2 << std::endl;//空行
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
	//file2 << std::endl;//空行
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

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
