// NetworkAdjust.h: interface for the CNetworkAdjustment class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_NETWORKADJUST_H__029E7110_362B_45E9_BCA6_4D307D91111B__INCLUDED_)
#define AFX_NETWORKADJUST_H__029E7110_362B_45E9_BCA6_4D307D91111B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <Eigen/Dense>
#include <conio.h>
#include <iostream>
#include <math.h>
#include <optional>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class CPlaneNetAdjust {
public:
    FILE* resultfp = nullptr; // 结果文件指针
    int m_Pnumber = 0; // 总点数
    int m_Lnumber = 0; // 方向组总数
    int m_Nnumber = 0; // 方向值总数
    int m_knPnumber = 0; // 已知点总数
    int m_StableNumber = 0; // 拟稳点数

    int m_Tnumber = 0; // 方位角观测值个数
    int m_Snumber = 0; // 边长观测值个数

    int* dir1 = nullptr; // 测站点号
    int* dir2 = nullptr; // 观测方向的点号
    int* dir0 = nullptr; // 测站零方向在方向值数组中的位置

    char** Pname = nullptr; // 点名地址数组

    double *L = nullptr, *V = nullptr; // 方向观测值及其残差数组
    double* XY = nullptr; // 坐标数组
    double *ATPA = nullptr, *ATPL = nullptr;
    // 法方程系数阵、法方程自由项
    Eigen::MatrixXd m_ATPA;
    Eigen::VectorXd m_ATPL;
    double *ATPA1 = nullptr, *ATPA2 = nullptr;
    // 角度分块法方程系数、距离分块法方程系数
    double m0a = 0.0, m0s = 0.0; // 初始测角中误差和测距中误差
    double* dX = nullptr; // 未知参数数组（坐标改正数、定向角改正数）
    double m_pvv = 0.0; // [pvv]
    double m_mu = 0.0; // 验后单位权中误差

    bool* IsStable = nullptr; // 拟稳点标志数组

    double ma = 0.0, mT = 0.0; // 方向值中误差、方位角中误差
    double mS1 = 0.0, mS2 = 0.0; // 边长固定误差、比例误差

    int* T_dir1 = nullptr; // 方位角测站点点号数组
    int* T_dir2 = nullptr; // 方位角照准点点号数组
    double* T_L = nullptr; // 方位角观测值数组
    double* T_V = nullptr; // 方位角残差数组

    int* S_dir1 = nullptr; // 边长观测值测站点号数组
    int* S_dir2 = nullptr; // 边长观测值照准点号数组
    double* S_L = nullptr; // 边长观测值数组
    double* S_V = nullptr; // 边长残差数组

    bool* Usable = nullptr; // 观测值可用标志数组

    double m_A_P = 0.0; // 方向观测值输入权
    double m_T_P = 0.0; // 方位角观测值输入权
    double m_S_P = 0.0; // 边长观测值输入权
    bool m_bA_P = false; // 是否输入方向观测值权
    bool m_bT_P = false; // 是否输入方位角观测值权
    bool m_bS_P = false; // 是否输入边长观测值权

    CPlaneNetAdjust();
    virtual ~CPlaneNetAdjust();

    void InputData(char* DATAFILE); // 输入原始数据
    void PrintData(); // 输出原始数据
    int GetStationNumber(char* name); // 保存点名，返回点号

    double Get_S12(int k1, int k2); // 查找边长观测值
    double Get_T12(int k1, int k2); // 查找方位角观测值
    double Get_Angle(int p, int p1, int p2); // 查找角度值
    void ca_xy0(); // 近似坐标计算(三角网)
    int cal_xy1(); // 近似坐标计算（自由平面网）
    int ca_x0y0(); // 近似坐标计算（导线网）

    double ca_T12(int k1, int k2); // 用坐标计算方位角，返回值为弧度值,
    double ca_ab(int k1, int k2, double B[], int Bin[]); // ab系数计算
    double ca_cd(int k1, int k2, double B[], int Bin[]); // cd系数计算
    void ca_ATPA(); // 法方程组成
    void ca_ATPA12(); // 法方程组成
    void ca_ATPAi(double B[], int Bin[], double p, double Li, int m); // 法方程累加项计算
    void ca_ATPA1(double B[], int Bin[], double p, double Li, int m); // 法方程累加项计算(Helmert中角度系数)
    void ca_ATPA2(double B[], int Bin[], double p, double Li, int m); // 法方程累加项计算（helmert中边长系数）
    double ca_dX(); // 参数平差值计算
    double ca_V(); // 残差计算
    double qq(double B[], int Bin[]); // 权倒数计算
    void PrintResult(double* mx, double* my, double* M,
        double* R_T, double* R_mT, double* R_S, double* R_ms,
        double* S_T, double* S_mT, double* S_SVs, double* S_ms,
        double* T_TV, double* T_mT, double* T_S, double* T_ms); // 输出平差成果

    int LS_Adjust(); // 最小二乘平差
    double LS_MaxValue; // 平面最小二乘平差改正数限值，默认0.01
    void Free(); // 自由网平差
    int Quasi_Stable(char* QuasiStablefile); // 拟稳平差
    double* ca_GT(char* file); // 自由网平差用的G矩阵
    int Helmert_LS(); // helmert方差分量估计
    int m_NHelmertCircle; // helmert方差分量估计循环次数，默认10次

    void ErrorEllipse(); // 误差椭圆计算

    double Snooping(double arfa);
};

#endif // !defined(AFX_NETWORKADJUST_H__029E7110_362B_45E9_BCA6_4D307D91111B__INCLUDED_)