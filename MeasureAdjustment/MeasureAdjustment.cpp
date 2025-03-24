#include "framework.h"
#include "pch.h"

#include "MeasureAdjustment.h"
// #include <StdAfx.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

//-----------水准网平差-----------
// 单位：高程：米，距离：千米
// 输入：
//		Lnumber:观测高差总数
//		Pnumber：总点数
//		knPnumber：已知点数
//		Sigma:验前单位权中误差,每千米高差中误差，以米为单位。
//
//		Pname:点名数组，所有点名存储在数组中，其在数组中的索引对应其实际序号
//		Height：已知高程值数组，其索引和相应的点名序号对应,该参数也是输出参数
//		StartP:观测数据  高差起点号（里面存储点名对应的序号）数组
//		EndP:观测数据  高差终点号（里面存储点名对应的序号）数组
//		L,P:观测数据  高差值L数组，路线长度P数组，该数组索引和StartP和EndP数组索引对应
// 返回值：
//      pvv:平差结果pvv
//      mu:验后单位权中误差
//
//      Height：高程平差值（数组），序号和点名相对应
//      dx:高程改正数（数组），序号和点名想对应；Height-dx即为近似高程
//      m0:高程中误差（数组），序号和点名相对应
//
//      V:高差改正数（数组），序号和观测数据序号对应，L+V是高差平差值
//      P:观测高差段权（数组）
//      m1:高差中误差（数组），序号和观测数据序号对应
//
int LevelAD(const int Lnumber, const int Pnumber, const int knPnumber, const double Sigma, char** Pname, double* Height, int* StartP, int* EndP, double* L, double* P, double& pvv, double& mu, double* dx, double* m0, double* V, double* m1)
{
    CLevelingAdjust lvad;
    lvad.m_Pnumber = Pnumber;
    lvad.m_Lnumber = Lnumber;
    lvad.m_knPnumber = knPnumber;
    lvad.m_Sigma = Sigma;
    lvad.Height = new double[Pnumber];
    lvad.Pname = new char*[Pnumber];
    lvad.dX = new double[lvad.m_Pnumber];
    lvad.ATPA = new double[lvad.m_Pnumber * (lvad.m_Pnumber + 1) / 2];
    lvad.ATPL = new double[lvad.m_Pnumber];
    lvad.StartP = new int[Lnumber];
    lvad.EndP = new int[Lnumber];
    lvad.L = new double[Lnumber];
    lvad.V = new double[Lnumber];
    lvad.P = new double[Lnumber];
    for (size_t i = 0; i < Pnumber; i++) {
        int len = strlen(Pname[i]);
        lvad.Pname[i] = new char[len + 1];
        strcpy_s(lvad.Pname[i], len + 1, Pname[i]);
        lvad.Height[i] = Height[i];
    }
    for (size_t i = 0; i < Lnumber; i++) {
        P[i] = 1.0 / P[i]; // 变为权
        lvad.StartP[i] = StartP[i];
        lvad.EndP[i] = EndP[i];
        lvad.L[i] = L[i];
        // lvad.V[i] = V[i];
        lvad.P[i] = P[i];
    }

    lvad.LS_Adjustment();
    pvv = lvad.m_pvv;
    mu = lvad.m_mu;
    for (size_t i = 0; i < Pnumber; i++) {
        Height[i] = lvad.Height[i];
        dx[i] = lvad.dX[i];
        double qii = lvad.ATPA[ij(i, i)];
        m0[i] = sqrt(qii) * mu;
    }
    for (size_t i = 0; i < Lnumber; i++) {
        int k1 = StartP[i];
        int k2 = EndP[i];
        double qii = lvad.ATPA[ij(k1, k1)];
        double qjj = lvad.ATPA[ij(k2, k2)];
        double qij = lvad.ATPA[ij(k1, k2)];
        m1[i] = sqrt(qii + qjj - 2.0 * qij) * mu;
        V[i] = lvad.V[i];
    }

    return 1;
}

// 水准网附合路线闭合差计算closure error
// 单位：高程：米，距离：千米
// 输入：
//		Lnumber:观测高差总数
//		Pnumber：总点数
//		knPnumber：已知点数
//		Sigma:验前单位权中误差
//
//		Pname:点名数组，所有点名存储在数组中，其在数组中的索引对应其实际序号
//		Height：已知高程值数组，其索引和相应的点名序号对应,该参数也是输出参数
//		StartP:观测数据  高差起点号（里面存储点名对应的序号）数组
//		EndP:观测数据  高差终点号（里面存储点名对应的序号）数组
//		L,P:观测数据  高差值L数组，路线长度P数组，该数组索引和StartP和EndP数组索引对应
// 返回值：
//       line_name:附合路线名称数组，比如line_name[1]为A B C D F，点名与点名之间为空格
//       line_L:符合路线长度数组
//       line_w:相应符合路线的闭合差数组
//       line_limit:相应符合路线的限差数组
//       loop_name:环闭合路线名称数组，如loop_name[1]为 A B C A,点名与点名之间为空格
//       loop_L:闭合环路线长度
//       loop_w:闭合差
//       loop_limit:限差
int LevelCE(const int Lnumber, const int Pnumber, const int knPnumber, const double Sigma, char** Pname, double* Height, int* StartP, int* EndP, double* L, double* P,
    std::vector<std::string>& line_name, std::vector<double>& line_L, std::vector<double>& line_w, std::vector<double>& line_limit,
    std::vector<std::string>& loop_name, std::vector<double>& loop_L, std::vector<double>& loop_w, std::vector<double>& loop_limit)
{
    CLevelingAdjust lvad;
    lvad.m_Pnumber = Pnumber;
    lvad.m_Lnumber = Lnumber;
    lvad.m_knPnumber = knPnumber;
    lvad.m_Sigma = Sigma;
    lvad.Height = new double[Pnumber];
    lvad.Pname = new char*[Pnumber];
    lvad.dX = new double[Pnumber];
    lvad.ATPA = new double[Pnumber * (Pnumber + 1) / 2];
    lvad.ATPL = new double[Pnumber];
    lvad.StartP = new int[Lnumber];
    lvad.EndP = new int[Lnumber];
    lvad.L = new double[Lnumber];
    lvad.V = new double[Lnumber];
    lvad.P = new double[Lnumber];
    for (size_t i = 0; i < Pnumber; i++) {
        int len = strlen(Pname[i]);
        lvad.Pname[i] = new char[len + 1];
        strcpy_s(lvad.Pname[i], len + 1, Pname[i]);
        lvad.Height[i] = Height[i];
    }
    for (size_t i = 0; i < Lnumber; i++) {
        P[i] = 1.0 / P[i]; // 变为权
        lvad.StartP[i] = StartP[i];
        lvad.EndP[i] = EndP[i];
        lvad.L[i] = L[i];
        lvad.P[i] = P[i];
    }

    lvad.LineClosure(line_name, line_L, line_w, line_limit);
    lvad.LoopClosure(loop_name, loop_L, loop_w, loop_limit);

    return 1;
}

// 水平网最小二乘平差
//  单位：角度和方向观测数据为弧度，边长观测单位为米
//  输入：
//       Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:总点数、已知点数、方向组总数、方向值总数、边长观测值个数、方位角观测值个数
//       ma,mS1,mS2,mT:方向值中误差(以秒为单位)、边长固定误差（以米为单位）、边长的比例误差（无量纲）、方位角中误差（以秒为单位）
//       Pname:点名数组，大小为所有点大小
//       dir0,dir1,dir2,L:存储方向值观测数据；dir0方向组数组，大小方向组总数+1，存储序号为每组方向数据第一行观测值在L中存储的序号；
//                       dir1测站点好，大小为方向组总数；dir2照准点号，大小为方向值总数；L存储方向值，大小方向值总数。
//       S_dir1,S_dir2,S_L：存储边长观测数据，大小均为边长观测值总数；S_dir1测站点号，S_dir2照准点号，S_L边长观测值。
//       T_dir1,T_dir2,T_L:存储方位角观测数据，大小均为方位角观测总数；T_dir1测站点号，T_dir2照准点号，T_L方位角观测值。
//       b_manuAP,b_manu_SP,b_manu_TP:是否输入方向值、边长、方位角权重
//       manuAP,manuSP,manuTP:输入的方向值、边长、方位角权重
//       LS_max:最小二乘平差中循环残差限差，单位是米
//
//  返回：V,S_V,T_V分别是方向值误差、边长误差和方位角误差
//       S_SVs为S+V是边长平差值；T_TV是T+V是方位角观测平差值
//       XY，mx,my,M是点坐标平差值及其误差
//       mu是网平差精度
//
int PlaneNetAD(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L, double* V, double* R_T, double* R_mT, double* R_S, double* R_ms,
    int* S_dir1, int* S_dir2, double* S_L, double* S_V, double* S_T, double* S_mT, double* S_SVs, double* S_ms,
    int* T_dir1, int* T_dir2, double* T_L, double* T_V, double* T_TV, double* T_mT, double* T_S, double* T_ms,
    double* XY, double* mx, double* my, double* M,
    double& mu, double LS_max, bool b_manuAP, double manuAP, bool b_manuSP, double manuSP, bool b_manuTP, double manuTP)
{
    CPlaneNetAdjust cpa;
    cpa.m_Pnumber = Pnumber;
    cpa.m_knPnumber = knPnumber;
    cpa.m_Lnumber = Lnumber;
    cpa.m_Nnumber = Nnumber;
    cpa.m_Snumber = Snumber;
    cpa.m_Tnumber = Tnumber;

    cpa.ma = ma;
    cpa.mS1 = mS1;
    cpa.mS2 = mS2;
    cpa.mT = mT;

    cpa.m_bA_P = b_manuAP;
    cpa.m_bS_P = b_manuSP;
    cpa.m_bT_P = b_manuTP;
    cpa.m_A_P = manuAP;
    cpa.m_S_P = manuSP;
    cpa.m_T_P = manuTP;

    cpa.XY = new double[2 * Pnumber];
    cpa.Pname = new char*[Pnumber];
    for (size_t i = 0; i < Pnumber; i++) // 点名
    {
        int len = strlen(Pname[i]);
        cpa.Pname[i] = new char[len + 1];
        strcpy_s(cpa.Pname[i], len + 1, Pname[i]);
    }
    if (Lnumber > 0) // 方向观测值
    {
        cpa.dir0 = new int[Lnumber + 1];
        cpa.dir1 = new int[Lnumber];
        cpa.dir2 = new int[Nnumber];
        cpa.L = new double[Nnumber];
        cpa.V = new double[Nnumber];
        for (size_t i = 0; i <= Lnumber; i++) {
            cpa.dir0[i] = dir0[i];
            if (i != Lnumber) {
                cpa.dir1[i] = dir1[i];
            }
        }
        for (size_t i = 0; i < Nnumber; i++) {
            cpa.dir2[i] = dir2[i];
            cpa.L[i] = L[i];
        }
    }
    if (Snumber > 0) // 边长观测值
    {
        cpa.S_dir1 = new int[Snumber]; // 测站点号
        cpa.S_dir2 = new int[Snumber]; // 照准点号
        cpa.S_L = new double[Snumber]; // 边长观测值
        cpa.S_V = new double[Snumber]; // 边长残差
        for (size_t i = 0; i < Snumber; i++) {
            cpa.S_dir1[i] = S_dir1[i];
            cpa.S_dir2[i] = S_dir2[i];
            cpa.S_L[i] = S_L[i];
        }
    }
    if (Tnumber > 0) // 方位角观测值
    {
        cpa.T_dir1 = new int[Tnumber]; // 测站点号
        cpa.T_dir2 = new int[Tnumber]; // 照准点号
        cpa.T_L = new double[Tnumber]; // 方位角观测值
        cpa.T_V = new double[Tnumber]; // 方位角残差
        for (size_t i = 0; i < Tnumber; i++) {
            cpa.T_dir1[i] = T_dir1[i];
            cpa.T_dir2[i] = T_dir2[i];
            cpa.T_L[i] = T_L[i];
        }
    }
    // 观测值可用标志数组
    int n = cpa.m_Nnumber + cpa.m_Snumber + cpa.m_Tnumber;
    cpa.Usable = new bool[n];
    for (int i = 0; i < n; i++)
        cpa.Usable[i] = true;

    int t = 2 * cpa.m_Pnumber + cpa.m_Lnumber; // 未知参数总数
    int tt = t * (t + 1) / 2;

    cpa.ATPL = new double[t]; // 法方程自由项
    cpa.ATPA = new double[tt]; // 系数矩阵
    cpa.dX = new double[t]; // 未知数向量

    int unPnumber = cpa.m_Pnumber - cpa.m_knPnumber;

    for (size_t i = 0; i < 2 * knPnumber; i++) {
        cpa.XY[i] = XY[i];
    }

    // cpa.ca_x0y0();//近视坐标计算
    cpa.LS_MaxValue = LS_max;
    if (cpa.ca_x0y0() != 1)
        return 0; // 近似坐标计算
    if (cpa.LS_Adjust() != 1)
        return 0; // 最小二乘计算
    for (size_t i = 0; i < 2 * Pnumber; i++) {
        XY[i] = cpa.XY[i];
    }
    for (size_t i = 0; i < Nnumber; i++) {
        V[i] = cpa.V[i];
    }
    for (size_t i = 0; i < Snumber; i++) {
        S_V[i] = cpa.S_V[i];
    }
    for (size_t i = 0; i < Tnumber; i++) {
        T_V[i] = cpa.T_V[i];
    }
    mu = cpa.m_mu;
    cpa.PrintResult(mx, my, M,
        R_T, R_mT, R_S, R_ms,
        S_T, S_mT, S_SVs, S_ms,
        T_TV, T_mT, T_S, T_ms); // 结果精度计算

    return 1;
}

// 水平网最小二乘平差（专用，初始点计算方式不同）
//  单位：角度和方向观测数据为弧度，边长观测单位为米
//  输入：
//       Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:总点数、已知点数、方向组总数、方向值总数、边长观测值个数、方位角观测值个数
//       ma,mS1,mS2,mT:方向值中误差(以秒为单位)、边长固定误差（以米为单位）、边长的比例误差（无量纲）、方位角中误差（以秒为单位）
//       Pname:点名数组，大小为所有点大小
//       dir0,dir1,dir2,L:存储方向值观测数据；dir0方向组数组，大小方向组总数+1，存储序号为每组方向数据第一行观测值在L中存储的序号；
//                       dir1测站点好，大小为方向组总数；dir2照准点号，大小为方向值总数；L存储方向值，大小方向值总数。
//       S_dir1,S_dir2,S_L：存储边长观测数据，大小均为边长观测值总数；S_dir1测站点号，S_dir2照准点号，S_L边长观测值。
//       T_dir1,T_dir2,T_L:存储方位角观测数据，大小均为方位角观测总数；T_dir1测站点号，T_dir2照准点号，T_L方位角观测值。
//
//  返回：V,S_V,T_V分别是方向值误差、边长误差和方位角误差
//       S_SVs为S+V是边长平差值；T_TV是T+V是方位角观测平差值
//       XY，mx,my,M是点坐标平差值及其误差
//       mu是网平差精度
//
int PlaneNetAD1(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L, double* V, double* R_T, double* R_mT, double* R_S, double* R_ms,
    int* S_dir1, int* S_dir2, double* S_L, double* S_V, double* S_T, double* S_mT, double* S_SVs, double* S_ms,
    int* T_dir1, int* T_dir2, double* T_L, double* T_V, double* T_TV, double* T_mT, double* T_S, double* T_ms,
    double* XY, double* mx, double* my, double* M,
    double& mu)
{
    CPlaneNetAdjust cpa;
    cpa.m_Pnumber = Pnumber;
    cpa.m_knPnumber = knPnumber;
    cpa.m_Lnumber = Lnumber;
    cpa.m_Nnumber = Nnumber;
    cpa.m_Snumber = Snumber;
    cpa.m_Tnumber = Tnumber;

    cpa.ma = ma;
    cpa.mS1 = mS1;
    cpa.mS2 = mS2;
    cpa.mT = mT;

    cpa.XY = new double[2 * Pnumber];
    cpa.Pname = new char*[Pnumber];
    for (size_t i = 0; i < Pnumber; i++) // 点名
    {
        int len = strlen(Pname[i]);
        cpa.Pname[i] = new char[len + 1];
        strcpy_s(cpa.Pname[i], len + 1, Pname[i]);
    }
    if (Lnumber > 0) // 方向观测值
    {
        cpa.dir0 = new int[Lnumber + 1];
        cpa.dir1 = new int[Lnumber];
        cpa.dir2 = new int[Nnumber];
        cpa.L = new double[Nnumber];
        cpa.V = new double[Nnumber];
        for (size_t i = 0; i <= Lnumber; i++) {
            cpa.dir0[i] = dir0[i];
            if (i != Lnumber) {
                cpa.dir1[i] = dir1[i];
            }
        }
        for (size_t i = 0; i < Nnumber; i++) {
            cpa.dir2[i] = dir2[i];
            cpa.L[i] = L[i];
        }
    }
    if (Snumber > 0) // 边长观测值
    {
        cpa.S_dir1 = new int[Snumber]; // 测站点号
        cpa.S_dir2 = new int[Snumber]; // 照准点号
        cpa.S_L = new double[Snumber]; // 边长观测值
        cpa.S_V = new double[Snumber]; // 边长残差
        for (size_t i = 0; i < Snumber; i++) {
            cpa.S_dir1[i] = S_dir1[i];
            cpa.S_dir2[i] = S_dir2[i];
            cpa.S_L[i] = S_L[i];
        }
    }
    if (Tnumber > 0) // 方位角观测值
    {
        cpa.T_dir1 = new int[Tnumber]; // 测站点号
        cpa.T_dir2 = new int[Tnumber]; // 照准点号
        cpa.T_L = new double[Tnumber]; // 方位角观测值
        cpa.T_V = new double[Tnumber]; // 方位角残差
        for (size_t i = 0; i < Tnumber; i++) {
            cpa.T_dir1[i] = T_dir1[i];
            cpa.T_dir2[i] = T_dir2[i];
            cpa.T_L[i] = T_L[i];
        }
    }
    // 观测值可用标志数组
    int n = cpa.m_Nnumber + cpa.m_Snumber + cpa.m_Tnumber;
    cpa.Usable = new bool[n];
    for (int i = 0; i < n; i++)
        cpa.Usable[i] = true;

    int t = 2 * cpa.m_Pnumber + cpa.m_Lnumber; // 未知参数总数
    int tt = t * (t + 1) / 2;

    cpa.ATPL = new double[t]; // 法方程自由项
    cpa.ATPA = new double[tt]; // 系数矩阵
    cpa.dX = new double[t]; // 未知数向量

    int unPnumber = cpa.m_Pnumber - cpa.m_knPnumber;

    // cpa.ca_x0y0();//近视坐标计算
    // 近视坐标计算
    if (cpa.cal_xy1()) {

        // 将计算的近视坐标输出
        std::ofstream file2("..\\计算的点坐标.txt");
        // 检查文件是否成功打开
        if (!file2.is_open()) {
            std::cerr << "无法打开文件 " << std::endl;
            return 0; // 如果文件打开失败，退出程序
        }
        for (size_t i = 0; i < Pnumber; i++) {
            file2 << Pname[i] << "\t" << cpa.XY[2 * i] << "\t" << cpa.XY[2 * i + 1] << std::endl;
        }
        file2.close();

        // cpa.m_knPnumber = cpa.m_Pnumber;
        // cpa.Free();//自由网平差

        double temp_dx, temp_dy, temp_s, temp_ro;
        calculateTransformationParameters(knPnumber, cpa.XY, XY, temp_dx, temp_dy, temp_s, temp_ro);
        // 将坐标转换到已知坐标的坐标系下
        for (size_t i = knPnumber; i < Pnumber; i++) {
            transformPoint(cpa.XY[2 * i], cpa.XY[2 * i + 1], XY[2 * i], XY[2 * i + 1], temp_dx, temp_dy, temp_s, temp_ro);
        }
    }
    for (size_t i = 0; i < 2 * knPnumber; i++) // 然后将所有的XY给复制过去
    {
        cpa.XY[i] = XY[i];
    }

    if (cpa.LS_Adjust() != 1)
        return 0; // 最小二乘计算
    for (size_t i = 0; i < 2 * Pnumber; i++) {
        XY[i] = cpa.XY[i];
    }
    for (size_t i = 0; i < Nnumber; i++) {
        V[i] = cpa.V[i];
    }
    for (size_t i = 0; i < Snumber; i++) {
        S_V[i] = cpa.S_V[i];
    }
    for (size_t i = 0; i < Tnumber; i++) {
        T_V[i] = cpa.T_V[i];
    }
    mu = cpa.m_mu;
    cpa.PrintResult(mx, my, M,
        R_T, R_mT, R_S, R_ms,
        S_T, S_mT, S_SVs, S_ms,
        T_TV, T_mT, T_S, T_ms); // 结果精度计算

    return 1;
}

// 水平网自由网平差
//  单位：角度和方向观测数据为弧度，边长观测单位为米
//  输入：
//       Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:总点数、已知点数、方向组总数、方向值总数、边长观测值个数、方位角观测值个数
//       ma,mS1,mS2,mT:方向值中误差、方位角中误差、边长固定误差、比例误差
//       Pname:点名数组，大小为所有点大小
//       dir0,dir1,dir2,L:存储方向值观测数据；dir0方向组数组，大小方向组总数+1，存储序号为每组方向数据第一行观测值在L中存储的序号；
//                       dir1测站点好，大小为方向组总数；dir2照准点号，大小为方向值总数；L存储方向值，大小方向值总数。
//       S_dir1,S_dir2,S_L：存储边长观测数据，大小均为边长观测值总数；S_dir1测站点号，S_dir2照准点号，S_L边长观测值。
//       T_dir1,T_dir2,T_L:存储方位角观测数据，大小均为方位角观测总数；T_dir1测站点号，T_dir2照准点号，T_L方位角观测值。
//
//  返回：V,S_V,T_V分别是方向值误差、边长误差和方位角误差
//       S_SVs为S+V是边长平差值；T_TV是T+V是方位角观测平差值
//       XY，mx,my,M是点坐标平差值及其误差
//       mu是网平差精度
//
int PlaneNetFree(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L, double* V, double* R_T, double* R_mT, double* R_S, double* R_ms,
    int* S_dir1, int* S_dir2, double* S_L, double* S_V, double* S_T, double* S_mT, double* S_SVs, double* S_ms,
    int* T_dir1, int* T_dir2, double* T_L, double* T_V, double* T_TV, double* T_mT, double* T_S, double* T_ms,
    double* XY, double* mx, double* my, double* M,
    double& mu)
{
    CPlaneNetAdjust cpfree;
    cpfree.m_Pnumber = Pnumber;
    cpfree.m_knPnumber = knPnumber;
    cpfree.m_Lnumber = Lnumber;
    cpfree.m_Nnumber = Nnumber;
    cpfree.m_Snumber = Snumber;
    cpfree.m_Tnumber = Tnumber;

    cpfree.ma = ma;
    cpfree.mS1 = mS1;
    cpfree.mS2 = mS2;
    cpfree.mT = mT;

    cpfree.XY = new double[2 * Pnumber];
    cpfree.Pname = new char*[Pnumber];
    for (size_t i = 0; i < Pnumber; i++) // 点名
    {
        int len = strlen(Pname[i]);
        cpfree.Pname[i] = new char[len + 1];
        strcpy_s(cpfree.Pname[i], len + 1, Pname[i]);
    }
    if (Lnumber > 0) // 方向观测值
    {
        cpfree.dir0 = new int[Lnumber + 1];
        cpfree.dir1 = new int[Lnumber];
        cpfree.dir2 = new int[Nnumber];
        cpfree.L = new double[Nnumber];
        cpfree.V = new double[Nnumber];
        for (size_t i = 0; i <= Lnumber; i++) {
            cpfree.dir0[i] = dir0[i];
            if (i != Lnumber) {
                cpfree.dir1[i] = dir1[i];
            }
        }
        for (size_t i = 0; i < Nnumber; i++) {
            cpfree.dir2[i] = dir2[i];
            cpfree.L[i] = L[i];
        }
    }
    if (Snumber > 0) // 边长观测值
    {
        cpfree.S_dir1 = new int[Snumber]; // 测站点号
        cpfree.S_dir2 = new int[Snumber]; // 照准点号
        cpfree.S_L = new double[Snumber]; // 边长观测值
        cpfree.S_V = new double[Snumber]; // 边长残差
        for (size_t i = 0; i < Snumber; i++) {
            cpfree.S_dir1[i] = S_dir1[i];
            cpfree.S_dir2[i] = S_dir2[i];
            cpfree.S_L[i] = S_L[i];
        }
    }
    if (Tnumber > 0) // 方位角观测值
    {
        cpfree.T_dir1 = new int[Tnumber]; // 测站点号
        cpfree.T_dir2 = new int[Tnumber]; // 照准点号
        cpfree.T_L = new double[Tnumber]; // 方位角观测值
        cpfree.T_V = new double[Tnumber]; // 方位角残差
        for (size_t i = 0; i < Tnumber; i++) {
            cpfree.T_dir1[i] = T_dir1[i];
            cpfree.T_dir2[i] = T_dir2[i];
            cpfree.T_L[i] = T_L[i];
        }
    }
    // 观测值可用标志数组
    int n = cpfree.m_Nnumber + cpfree.m_Snumber + cpfree.m_Tnumber;
    cpfree.Usable = new bool[n];
    for (int i = 0; i < n; i++)
        cpfree.Usable[i] = true;

    int t = 2 * cpfree.m_Pnumber + cpfree.m_Lnumber; // 未知参数总数
    int tt = t * (t + 1) / 2;

    cpfree.ATPL = new double[t]; // 法方程自由项
    cpfree.ATPA = new double[tt]; // 系数矩阵
    cpfree.dX = new double[t]; // 未知数向量

    int unPnumber = cpfree.m_Pnumber - cpfree.m_knPnumber;

    // cpfree.ca_x0y0();//近视坐标计算
    cpfree.cal_xy1(); // 近视坐标计算
    for (size_t i = 0; i < Pnumber; i++) {
        std::cout << i << "\t" << cpfree.Pname[i] << "\t" << cpfree.XY[2 * i] << "\t" << cpfree.XY[2 * i + 1] << std::endl;
    }

    cpfree.Free(); // 最小二乘计算
    for (size_t i = 0; i < 2 * Pnumber; i++) {
        XY[i] = cpfree.XY[i];
    }
    for (size_t i = 0; i < Nnumber; i++) {
        V[i] = cpfree.V[i];
    }
    for (size_t i = 0; i < Snumber; i++) {
        S_V[i] = cpfree.S_V[i];
    }
    for (size_t i = 0; i < Tnumber; i++) {
        T_V[i] = cpfree.T_V[i];
    }
    mu = cpfree.m_mu;
    cpfree.PrintResult(mx, my, M,
        R_T, R_mT, R_S, R_ms,
        S_T, S_mT, S_SVs, S_ms,
        T_TV, T_mT, T_S, T_ms); // 结果精度计算

    return 1;
}

// 水平网Helmert方差分量估计（最小二乘平差）
//  单位：角度和方向观测数据为弧度，边长观测单位为米
//  输入：
//       Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:总点数、已知点数、方向组总数、方向值总数、边长观测值个数、方位角观测值个数
//       ma,mS1,mS2,mT:方向值中误差(以秒为单位)、边长固定误差（以米为单位）、边长的比例误差（无量纲）、方位角中误差（以秒为单位）
//       Pname:点名数组，大小为所有点大小
//       dir0,dir1,dir2,L:存储方向值观测数据；dir0方向组数组，大小方向组总数+1，存储序号为每组方向数据第一行观测值在L中存储的序号；
//                       dir1测站点好，大小为方向组总数；dir2照准点号，大小为方向值总数；L存储方向值，大小方向值总数。
//       S_dir1,S_dir2,S_L：存储边长观测数据，大小均为边长观测值总数；S_dir1测站点号，S_dir2照准点号，S_L边长观测值。
//       T_dir1,T_dir2,T_L:存储方位角观测数据，大小均为方位角观测总数；T_dir1测站点号，T_dir2照准点号，T_L方位角观测值。
//
//  返回：V,S_V,T_V分别是方向值误差、边长误差和方位角误差
//       S_SVs为S+V是边长平差值；T_TV是T+V是方位角观测平差值
//       XY，mx,my,M是点坐标平差值及其误差
//       mu是网平差精度,m1是helmert方差分量估计后的角度误差，m2是helmert方差分量估计后的距离误差
//       HelCircle:Helmert循环最大次数
//
int PlaneNetHelmertAd(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L, double* V, double* R_T, double* R_mT, double* R_S, double* R_ms,
    int* S_dir1, int* S_dir2, double* S_L, double* S_V, double* S_T, double* S_mT, double* S_SVs, double* S_ms,
    int* T_dir1, int* T_dir2, double* T_L, double* T_V, double* T_TV, double* T_mT, double* T_S, double* T_ms,
    double* XY, double* mx, double* my, double* M,
    double& mu, double& m1, double& m2, int HelCircle)
{
    CPlaneNetAdjust cpa;
    cpa.m_Pnumber = Pnumber;
    cpa.m_knPnumber = knPnumber;
    cpa.m_Lnumber = Lnumber;
    cpa.m_Nnumber = Nnumber;
    cpa.m_Snumber = Snumber;
    cpa.m_Tnumber = Tnumber;

    cpa.ma = ma;
    cpa.mS1 = mS1;
    cpa.mS2 = mS2;
    cpa.mT = mT;

    cpa.XY = new double[2 * Pnumber];
    cpa.Pname = new char*[Pnumber];
    for (size_t i = 0; i < Pnumber; i++) // 点名
    {
        int len = strlen(Pname[i]);
        cpa.Pname[i] = new char[len + 1];
        strcpy_s(cpa.Pname[i], len + 1, Pname[i]);
    }
    if (Lnumber > 0) // 方向观测值
    {
        cpa.dir0 = new int[Lnumber + 1];
        cpa.dir1 = new int[Lnumber];
        cpa.dir2 = new int[Nnumber];
        cpa.L = new double[Nnumber];
        cpa.V = new double[Nnumber];
        for (size_t i = 0; i <= Lnumber; i++) {
            cpa.dir0[i] = dir0[i];
            if (i != Lnumber) {
                cpa.dir1[i] = dir1[i];
            }
        }
        for (size_t i = 0; i < Nnumber; i++) {
            cpa.dir2[i] = dir2[i];
            cpa.L[i] = L[i];
        }
    }
    if (Snumber > 0) // 边长观测值
    {
        cpa.S_dir1 = new int[Snumber]; // 测站点号
        cpa.S_dir2 = new int[Snumber]; // 照准点号
        cpa.S_L = new double[Snumber]; // 边长观测值
        cpa.S_V = new double[Snumber]; // 边长残差
        for (size_t i = 0; i < Snumber; i++) {
            cpa.S_dir1[i] = S_dir1[i];
            cpa.S_dir2[i] = S_dir2[i];
            cpa.S_L[i] = S_L[i];
        }
    }
    if (Tnumber > 0) // 方位角观测值
    {
        cpa.T_dir1 = new int[Tnumber]; // 测站点号
        cpa.T_dir2 = new int[Tnumber]; // 照准点号
        cpa.T_L = new double[Tnumber]; // 方位角观测值
        cpa.T_V = new double[Tnumber]; // 方位角残差
        for (size_t i = 0; i < Tnumber; i++) {
            cpa.T_dir1[i] = T_dir1[i];
            cpa.T_dir2[i] = T_dir2[i];
            cpa.T_L[i] = T_L[i];
        }
    }
    // 观测值可用标志数组
    int n = cpa.m_Nnumber + cpa.m_Snumber + cpa.m_Tnumber;
    cpa.Usable = new bool[n];
    for (int i = 0; i < n; i++)
        cpa.Usable[i] = true;

    int t = 2 * cpa.m_Pnumber + cpa.m_Lnumber; // 未知参数总数
    int tt = t * (t + 1) / 2;

    cpa.ATPL = new double[t]; // 法方程自由项
    cpa.ATPA = new double[tt]; // 系数矩阵
    cpa.ATPA1 = new double[tt];
    cpa.ATPA2 = new double[tt];
    cpa.dX = new double[t]; // 未知数向量

    int unPnumber = cpa.m_Pnumber - cpa.m_knPnumber;

    for (size_t i = 0; i < 2 * knPnumber; i++) {
        cpa.XY[i] = XY[i];
    }

    // cpa.ca_x0y0();//近视坐标计算
    cpa.m_NHelmertCircle = HelCircle;
    if (cpa.ca_x0y0() != 1)
        return 0; // 近似坐标计算
    if (cpa.Helmert_LS() != 1)
        return 0; // Helmert方差分量估计
    m1 = cpa.m0a;
    m2 = cpa.m0s;
    for (size_t i = 0; i < 2 * Pnumber; i++) {
        XY[i] = cpa.XY[i];
    }
    for (size_t i = 0; i < Nnumber; i++) {
        V[i] = cpa.V[i];
    }
    for (size_t i = 0; i < Snumber; i++) {
        S_V[i] = cpa.S_V[i];
    }
    for (size_t i = 0; i < Tnumber; i++) {
        T_V[i] = cpa.T_V[i];
    }
    mu = cpa.m_mu;
    cpa.PrintResult(mx, my, M,
        R_T, R_mT, R_S, R_ms,
        S_T, S_mT, S_SVs, S_ms,
        T_TV, T_mT, T_S, T_ms); // 结果精度计算

    return 1;
}

// 水平网自由网初值计算
//  单位：角度和方向观测数据为弧度，边长观测单位为米
//  输入：
//       Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:总点数、已知点数、方向组总数、方向值总数、边长观测值个数、方位角观测值个数
//       ma,mS1,mS2,mT:方向值中误差(以秒为单位)、边长固定误差（以米为单位）、边长的比例误差（无量纲）、方位角中误差（以秒为单位）
//       Pname:点名数组，大小为所有点大小
//       dir0,dir1,dir2,L:存储方向值观测数据；dir0方向组数组，大小方向组总数+1，存储序号为每组方向数据第一行观测值在L中存储的序号；
//                       dir1测站点好，大小为方向组总数；dir2照准点号，大小为方向值总数；L存储方向值，大小方向值总数。
//       S_dir1,S_dir2,S_L：存储边长观测数据，大小均为边长观测值总数；S_dir1测站点号，S_dir2照准点号，S_L边长观测值。
//       T_dir1,T_dir2,T_L:存储方位角观测数据，大小均为方位角观测总数；T_dir1测站点号，T_dir2照准点号，T_L方位角观测值。
//
//  返回：
//       XY:计算的初始点坐标
//
int PlaneNetFreeinitPoint(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L,
    int* S_dir1, int* S_dir2, double* S_L,
    int* T_dir1, int* T_dir2, double* T_L,
    double* XY)
{
    CPlaneNetAdjust cpa;
    cpa.m_Pnumber = Pnumber;
    cpa.m_knPnumber = knPnumber;
    cpa.m_Lnumber = Lnumber;
    cpa.m_Nnumber = Nnumber;
    cpa.m_Snumber = Snumber;
    cpa.m_Tnumber = Tnumber;

    cpa.ma = ma;
    cpa.mS1 = mS1;
    cpa.mS2 = mS2;
    cpa.mT = mT;

    cpa.XY = new double[2 * Pnumber];
    cpa.Pname = new char*[Pnumber];
    for (size_t i = 0; i < Pnumber; i++) // 点名
    {
        int len = strlen(Pname[i]);
        cpa.Pname[i] = new char[len + 1];
        strcpy_s(cpa.Pname[i], len + 1, Pname[i]);
    }
    if (Lnumber > 0) // 方向观测值
    {
        cpa.dir0 = new int[Lnumber + 1];
        cpa.dir1 = new int[Lnumber];
        cpa.dir2 = new int[Nnumber];
        cpa.L = new double[Nnumber];
        cpa.V = new double[Nnumber];
        for (size_t i = 0; i <= Lnumber; i++) {
            cpa.dir0[i] = dir0[i];
            if (i != Lnumber) {
                cpa.dir1[i] = dir1[i];
            }
        }
        for (size_t i = 0; i < Nnumber; i++) {
            cpa.dir2[i] = dir2[i];
            cpa.L[i] = L[i];
        }
    }
    if (Snumber > 0) // 边长观测值
    {
        cpa.S_dir1 = new int[Snumber]; // 测站点号
        cpa.S_dir2 = new int[Snumber]; // 照准点号
        cpa.S_L = new double[Snumber]; // 边长观测值
        cpa.S_V = new double[Snumber]; // 边长残差
        for (size_t i = 0; i < Snumber; i++) {
            cpa.S_dir1[i] = S_dir1[i];
            cpa.S_dir2[i] = S_dir2[i];
            cpa.S_L[i] = S_L[i];
        }
    }
    if (Tnumber > 0) // 方位角观测值
    {
        cpa.T_dir1 = new int[Tnumber]; // 测站点号
        cpa.T_dir2 = new int[Tnumber]; // 照准点号
        cpa.T_L = new double[Tnumber]; // 方位角观测值
        cpa.T_V = new double[Tnumber]; // 方位角残差
        for (size_t i = 0; i < Tnumber; i++) {
            cpa.T_dir1[i] = T_dir1[i];
            cpa.T_dir2[i] = T_dir2[i];
            cpa.T_L[i] = T_L[i];
        }
    }
    // 观测值可用标志数组
    int n = cpa.m_Nnumber + cpa.m_Snumber + cpa.m_Tnumber;
    cpa.Usable = new bool[n];
    for (int i = 0; i < n; i++)
        cpa.Usable[i] = true;

    int t = 2 * cpa.m_Pnumber + cpa.m_Lnumber; // 未知参数总数
    int tt = t * (t + 1) / 2;

    cpa.ATPL = new double[t]; // 法方程自由项
    cpa.ATPA = new double[tt]; // 系数矩阵
    cpa.ATPA1 = new double[tt];
    cpa.ATPA2 = new double[tt];
    cpa.dX = new double[t]; // 未知数向量

    int unPnumber = cpa.m_Pnumber - cpa.m_knPnumber;

    for (size_t i = 0; i < 2 * knPnumber; i++) {
        cpa.XY[i] = XY[i];
    }

    cpa.cal_xy1(); // 近视坐标计算

    for (size_t i = 0; i < 2 * Pnumber; i++) {
        XY[i] = cpa.XY[i];
    }

    return 1;
}

// 水平网初值计算
//  单位：角度和方向观测数据为弧度，边长观测单位为米
//  输入：
//       Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:总点数、已知点数、方向组总数、方向值总数、边长观测值个数、方位角观测值个数
//       ma,mS1,mS2,mT:方向值中误差(以秒为单位)、边长固定误差（以米为单位）、边长的比例误差（无量纲）、方位角中误差（以秒为单位）
//       Pname:点名数组，大小为所有点大小
//       dir0,dir1,dir2,L:存储方向值观测数据；dir0方向组数组，大小方向组总数+1，存储序号为每组方向数据第一行观测值在L中存储的序号；
//                       dir1测站点好，大小为方向组总数；dir2照准点号，大小为方向值总数；L存储方向值，大小方向值总数。
//       S_dir1,S_dir2,S_L：存储边长观测数据，大小均为边长观测值总数；S_dir1测站点号，S_dir2照准点号，S_L边长观测值。
//       T_dir1,T_dir2,T_L:存储方位角观测数据，大小均为方位角观测总数；T_dir1测站点号，T_dir2照准点号，T_L方位角观测值。
//
//  返回：
//       XY:计算的初始点坐标
//
int PlaneNetinitPoint(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L,
    int* S_dir1, int* S_dir2, double* S_L,
    int* T_dir1, int* T_dir2, double* T_L,
    double* XY)
{
    CPlaneNetAdjust cpa;
    cpa.m_Pnumber = Pnumber;
    cpa.m_knPnumber = knPnumber;
    cpa.m_Lnumber = Lnumber;
    cpa.m_Nnumber = Nnumber;
    cpa.m_Snumber = Snumber;
    cpa.m_Tnumber = Tnumber;

    cpa.ma = ma;
    cpa.mS1 = mS1;
    cpa.mS2 = mS2;
    cpa.mT = mT;

    cpa.XY = new double[2 * Pnumber];
    cpa.Pname = new char*[Pnumber];
    for (size_t i = 0; i < Pnumber; i++) // 点名
    {
        int len = strlen(Pname[i]);
        cpa.Pname[i] = new char[len + 1];
        strcpy_s(cpa.Pname[i], len + 1, Pname[i]);
    }
    if (Lnumber > 0) // 方向观测值
    {
        cpa.dir0 = new int[Lnumber + 1];
        cpa.dir1 = new int[Lnumber];
        cpa.dir2 = new int[Nnumber];
        cpa.L = new double[Nnumber];
        cpa.V = new double[Nnumber];
        for (size_t i = 0; i <= Lnumber; i++) {
            cpa.dir0[i] = dir0[i];
            if (i != Lnumber) {
                cpa.dir1[i] = dir1[i];
            }
        }
        for (size_t i = 0; i < Nnumber; i++) {
            cpa.dir2[i] = dir2[i];
            cpa.L[i] = L[i];
        }
    }
    if (Snumber > 0) // 边长观测值
    {
        cpa.S_dir1 = new int[Snumber]; // 测站点号
        cpa.S_dir2 = new int[Snumber]; // 照准点号
        cpa.S_L = new double[Snumber]; // 边长观测值
        cpa.S_V = new double[Snumber]; // 边长残差
        for (size_t i = 0; i < Snumber; i++) {
            cpa.S_dir1[i] = S_dir1[i];
            cpa.S_dir2[i] = S_dir2[i];
            cpa.S_L[i] = S_L[i];
        }
    }
    if (Tnumber > 0) // 方位角观测值
    {
        cpa.T_dir1 = new int[Tnumber]; // 测站点号
        cpa.T_dir2 = new int[Tnumber]; // 照准点号
        cpa.T_L = new double[Tnumber]; // 方位角观测值
        cpa.T_V = new double[Tnumber]; // 方位角残差
        for (size_t i = 0; i < Tnumber; i++) {
            cpa.T_dir1[i] = T_dir1[i];
            cpa.T_dir2[i] = T_dir2[i];
            cpa.T_L[i] = T_L[i];
        }
    }
    // 观测值可用标志数组
    int n = cpa.m_Nnumber + cpa.m_Snumber + cpa.m_Tnumber;
    cpa.Usable = new bool[n];
    for (int i = 0; i < n; i++)
        cpa.Usable[i] = true;

    int t = 2 * cpa.m_Pnumber + cpa.m_Lnumber; // 未知参数总数
    int tt = t * (t + 1) / 2;

    cpa.ATPL = new double[t]; // 法方程自由项
    cpa.ATPA = new double[tt]; // 系数矩阵
    cpa.ATPA1 = new double[tt];
    cpa.ATPA2 = new double[tt];
    cpa.dX = new double[t]; // 未知数向量

    int unPnumber = cpa.m_Pnumber - cpa.m_knPnumber;

    for (size_t i = 0; i < 2 * knPnumber; i++) {
        cpa.XY[i] = XY[i];
    }

    if (cpa.ca_x0y0() != 1)
        return 0; // 近似坐标计算

    for (size_t i = 0; i < 2 * Pnumber; i++) {
        XY[i] = cpa.XY[i];
    }

    return 1;
}

// 水平网环闭合差计算

std::vector<std::pair<std::vector<std::string>, float>> PlaneCE(const int Lnumber, const int Nnumber, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L)
{
    std::vector<std::pair<std::vector<std::string>, float>> result;

    if (Lnumber > 0) // 方向观测值
    {
        Graph g;
        for (size_t i = 0; i < Lnumber; i++) {
            for (size_t j = dir0[i]; j < dir0[i + 1]; j++) {
                g.addEdge(std::to_string(dir1[i]), std::to_string(dir2[j]));
            }
        }
        vector<vector<string>> loops = g.findMinimalCycles();

        for (size_t i = 0; i < loops.size(); i++) {
            vector<string> loop = loops[i];
            loop.pop_back(); // 删除最后一个元素，这里最后一个元素和第一个元素相同，如果不同注释该行
            double plusAngle = 0;
            for (size_t j = 0; j < loop.size(); j++) {
                string frontP, backP; // 前一个、后一个点
                double frontA = 0, backA = 0; // 往前、往后的方位角
                string curP = loop[j];
                if (j == 0) {
                    frontP = loop[loop.size() - 1];
                } else {
                    frontP = loop[j - 1];
                }
                if (j == loop.size() - 1) {
                    backP = loop[0];
                } else {
                    backP = loop[j + 1];
                }
                // 查找当前点
                for (size_t k = 0; k < Lnumber; k++) {
                    if (dir1[k] == stoi(curP)) {
                        for (size_t k1 = dir0[k]; k1 < dir0[k + 1]; k1++) {
                            if (std::stoi(frontP) == dir2[k1]) {
                                frontA = L[k1];
                            }
                            if (std::stoi(backP) == dir2[k1]) {
                                backA = L[k1];
                            }
                        }
                        break;
                    }
                }
                // 计算当前角度
                double tempA = (backA - frontA) > 0 ? (backA - frontA) : (backA - frontA + 2 * EIGEN_PI);
                plusAngle += tempA;
            }
            // 多边形内角和(n-2)*180
            double Aangleerror = plusAngle - (loop.size() - 2) * EIGEN_PI;
            if (Aangleerror > EIGEN_PI / 4.0) // plusAngle非内角和
            {
                Aangleerror = plusAngle - (loop.size() + 2) * EIGEN_PI;
            }

            result.push_back(make_pair(loop, Aangleerror));
        }
    }
    return result;
}