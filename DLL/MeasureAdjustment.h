/*如果定义了宏就设置了__declspec(dllexport)修饰符，若未定义则设置__declspec(dllimport)*/
#ifdef MEASUREADJUSTMENT_EXPORTS
#define MEASUREADJUSTMENT_API __declspesc(dllexport)
#else
#define MEASUREADJUSTMENT_API __declspec(dllimport)
#endif // MEASUREADJUSTMENT


#include "pch.h"

#include "LevelingAdjust.h"
#include "public.h"
#include "PlaneNetAdjust.h"
#include <string>
#include <vector>
#include "Graph1.hpp"


//-----------水准网平差-----------
//单位：高程：米，距离：千米
//输入：
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
//返回值：
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
MEASUREADJUSTMENT_API int LevelAD(const int Lnumber, const int Pnumber, const int knPnumber, const double Sigma, char** Pname, double* Height, int* StartP, int* EndP, double* L, double* P, double& pvv, double& mu, double* dx, double* m0, double* V, double* m1);


//水准网附合路线闭合差计算closure error
//单位：高程：米，距离：千米
//输入：
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
//返回值：
//      line_name:附合路线名称数组，比如line_name[1]为A B C D F，点名与点名之间为空格
//      line_L:符合路线长度数组
//      line_w:相应符合路线的闭合差数组
//      line_limit:相应符合路线的限差数组
//      loop_name:环闭合路线名称数组，如loop_name[1]为 A B C A,点名与点名之间为空格
//      loop_L:闭合环路线长度
//      loop_w:闭合差
//      loop_limit:限差
MEASUREADJUSTMENT_API int LevelCE(const int Lnumber, const int Pnumber, const int knPnumber, const double Sigma, char** Pname, double* Height, int* StartP, int* EndP, double* L, double* P,
	std::vector<std::string>& line_name, std::vector<double>& line_L, std::vector<double>& line_w, std::vector<double>& line_limit,
	std::vector<std::string>& loop_name, std::vector<double>& loop_L, std::vector<double>& loop_w, std::vector<double>& loop_limit);

//水平网最小二乘平差(坐标单位为米)
// 单位：角度和方向观测数据为弧度，边长观测单位为米
// 输入：
//      Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:总点数、已知点数、方向组总数、方向值总数、边长观测值个数、方位角观测值个数
//      ma,mS1,mS2,mT:方向值中误差(以秒为单位)、边长固定误差（以米为单位）、边长的比例误差（无量纲）、方位角中误差（以秒为单位）
//      Pname:点名数组，大小为所有点大小
//      dir0,dir1,dir2,L:存储方向值观测数据；dir0方向组数组，大小方向组总数+1，存储序号为每组方向数据第一行观测值在L中存储的序号；
//                      dir1测站点好，大小为方向组总数；dir2照准点号，大小为方向值总数；L存储方向值，大小方向值总数。
//      S_dir1,S_dir2,S_L：存储边长观测数据，大小均为边长观测值总数；S_dir1测站点号，S_dir2照准点号，S_L边长观测值。
//      T_dir1,T_dir2,T_L:存储方位角观测数据，大小均为方位角观测总数；T_dir1测站点号，T_dir2照准点号，T_L方位角观测值。
//      b_manuAP,b_manu_SP,b_manu_TP:是否输入方向值、边长、方位角权重
//      manuAP,manuSP,manuTP:输入的方向值、边长、方位角权重
// //   LS_max:最小二乘平差中循环残差限差，单位是米，默认0.01
// 
// 返回：V,S_V,T_V分别是方向值误差、边长误差和方位角误差
//      S_SVs为S+V是边长平差值；T_TV是T+V是方位角观测平差值
//      XY，mx,my,M是点坐标平差值及其误差
//      mu是网平差精度
//
MEASUREADJUSTMENT_API int PlaneNetAD(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L, double* V, double* R_T, double* R_mT, double* R_S, double* R_ms,
    int* S_dir1, int* S_dir2, double* S_L, double* S_V, double* S_T, double* S_mT, double* S_SVs, double* S_ms,
    int* T_dir1, int* T_dir2, double* T_L, double* T_V, double* T_TV, double* T_mT, double* T_S, double* T_ms,
    double* XY, double* mx, double* my, double* M,
    double& mu, double LS_max=0.01, bool b_manuAP=false, double manuAP=1, bool b_manuSP=false, double manuSP=1, bool b_manuTP=false, double manuTP=1);
MEASUREADJUSTMENT_API int PlaneNetAD1(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L, double* V, double* R_T, double* R_mT, double* R_S, double* R_ms,
    int* S_dir1, int* S_dir2, double* S_L, double* S_V, double* S_T, double* S_mT, double* S_SVs, double* S_ms,
    int* T_dir1, int* T_dir2, double* T_L, double* T_V, double* T_TV, double* T_mT, double* T_S, double* T_ms,
    double* XY, double* mx, double* my, double* M,
    double& mu);

//水平网Helmert方差分量估计（最小二乘平差）
// 单位：角度和方向观测数据为弧度，边长观测单位为米
// 输入：
//      Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:总点数、已知点数、方向组总数、方向值总数、边长观测值个数、方位角观测值个数
//      ma,mS1,mS2,mT:方向值中误差(以秒为单位)、边长固定误差（以米为单位）、边长的比例误差（无量纲）、方位角中误差（以秒为单位）
//      Pname:点名数组，大小为所有点大小
//      dir0,dir1,dir2,L:存储方向值观测数据；dir0方向组数组，大小方向组总数+1，存储序号为每组方向数据第一行观测值在L中存储的序号；
//                      dir1测站点好，大小为方向组总数；dir2照准点号，大小为方向值总数；L存储方向值，大小方向值总数。
//      S_dir1,S_dir2,S_L：存储边长观测数据，大小均为边长观测值总数；S_dir1测站点号，S_dir2照准点号，S_L边长观测值。
//      T_dir1,T_dir2,T_L:存储方位角观测数据，大小均为方位角观测总数；T_dir1测站点号，T_dir2照准点号，T_L方位角观测值。
// 
// 返回：V,S_V,T_V分别是方向值误差、边长误差和方位角误差
//      S_SVs为S+V是边长平差值；T_TV是T+V是方位角观测平差值
//      XY，mx,my,M是点坐标平差值及其误差
//      mu是网平差精度
//      HelCircle:Helmert循环最大次数,默认10次
//
MEASUREADJUSTMENT_API int PlaneNetHelmertAd(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L, double* V, double* R_T, double* R_mT, double* R_S, double* R_ms,
    int* S_dir1, int* S_dir2, double* S_L, double* S_V, double* S_T, double* S_mT, double* S_SVs, double* S_ms,
    int* T_dir1, int* T_dir2, double* T_L, double* T_V, double* T_TV, double* T_mT, double* T_S, double* T_ms,
    double* XY, double* mx, double* my, double* M,
    double& mu, double& m1, double& m2,int HelCircle=10);

//水平网自由网平差
// 单位：角度和方向观测数据为弧度，边长观测单位为米
// 输入：
//      Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:总点数、已知点数、方向组总数、方向值总数、边长观测值个数、方位角观测值个数
//      ma,mS1,mS2,mT:方向值中误差(以秒为单位)、边长固定误差（以米为单位）、边长的比例误差（无量纲）、方位角中误差（以秒为单位）
//      Pname:点名数组，大小为所有点大小
//      dir0,dir1,dir2,L:存储方向值观测数据；dir0方向组数组，大小方向组总数+1，存储序号为每组方向数据第一行观测值在L中存储的序号；
//                      dir1测站点好，大小为方向组总数；dir2照准点号，大小为方向值总数；L存储方向值，大小方向值总数。
//      S_dir1,S_dir2,S_L：存储边长观测数据，大小均为边长观测值总数；S_dir1测站点号，S_dir2照准点号，S_L边长观测值。
//      T_dir1,T_dir2,T_L:存储方位角观测数据，大小均为方位角观测总数；T_dir1测站点号，T_dir2照准点号，T_L方位角观测值。
// 
// 返回：V,S_V,T_V分别是方向值误差、边长误差和方位角误差
//      S_SVs为S+V是边长平差值；T_TV是T+V是方位角观测平差值
//      XY，mx,my,M是点坐标平差值及其误差
//      mu是网平差精度
//
MEASUREADJUSTMENT_API int PlaneNetFree(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L, double* V, double* R_T, double* R_mT, double* R_S, double* R_ms,
    int* S_dir1, int* S_dir2, double* S_L, double* S_V, double* S_T, double* S_mT, double* S_SVs, double* S_ms,
    int* T_dir1, int* T_dir2, double* T_L, double* T_V, double* T_TV, double* T_mT, double* T_S, double* T_ms,
    double* XY, double* mx, double* my, double* M,
    double& mu);


//水平网初值计算
// 单位：角度和方向观测数据为弧度，边长观测单位为米
// 输入：
//      Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:总点数、已知点数、方向组总数、方向值总数、边长观测值个数、方位角观测值个数
//      ma,mS1,mS2,mT:方向值中误差(以秒为单位)、边长固定误差（以米为单位）、边长的比例误差（无量纲）、方位角中误差（以秒为单位）
//      Pname:点名数组，大小为所有点大小
//      dir0,dir1,dir2,L:存储方向值观测数据；dir0方向组数组，大小方向组总数+1，存储序号为每组方向数据第一行观测值在L中存储的序号；
//                      dir1测站点好，大小为方向组总数；dir2照准点号，大小为方向值总数；L存储方向值，大小方向值总数。
//      S_dir1,S_dir2,S_L：存储边长观测数据，大小均为边长观测值总数；S_dir1测站点号，S_dir2照准点号，S_L边长观测值。
//      T_dir1,T_dir2,T_L:存储方位角观测数据，大小均为方位角观测总数；T_dir1测站点号，T_dir2照准点号，T_L方位角观测值。
// 
// 返回：
//      XY:计算的初始点坐标
//
MEASUREADJUSTMENT_API int PlaneNetinitPoint(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L,
    int* S_dir1, int* S_dir2, double* S_L,
    int* T_dir1, int* T_dir2, double* T_L,
    double* XY);



//水平网自由网初值计算
// 单位：角度和方向观测数据为弧度，边长观测单位为米
// 输入：
//      Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:总点数、已知点数、方向组总数、方向值总数、边长观测值个数、方位角观测值个数
//      ma,mS1,mS2,mT:方向值中误差(以秒为单位)、边长固定误差（以米为单位）、边长的比例误差（无量纲）、方位角中误差（以秒为单位）
//      Pname:点名数组，大小为所有点大小
//      dir0,dir1,dir2,L:存储方向值观测数据；dir0方向组数组，大小方向组总数+1，存储序号为每组方向数据第一行观测值在L中存储的序号；
//                      dir1测站点好，大小为方向组总数；dir2照准点号，大小为方向值总数；L存储方向值，大小方向值总数。
//      S_dir1,S_dir2,S_L：存储边长观测数据，大小均为边长观测值总数；S_dir1测站点号，S_dir2照准点号，S_L边长观测值。
//      T_dir1,T_dir2,T_L:存储方位角观测数据，大小均为方位角观测总数；T_dir1测站点号，T_dir2照准点号，T_L方位角观测值。
// 
// 返回：
//      XY:计算的初始点坐标
//
MEASUREADJUSTMENT_API int PlaneNetFreeinitPoint(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L,
    int* S_dir1, int* S_dir2, double* S_L,
    int* T_dir1, int* T_dir2, double* T_L,
    double* XY);


//水平网环闭合差计算
//参数同上相同部分，返回指是一个pair组，包含组成环的点和环闭合差，闭合差单位是弧度
MEASUREADJUSTMENT_API std::vector<std::pair<std::vector<std::string>, float>> PlaneCE(const int Lnumber, const int Nnumber, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L);