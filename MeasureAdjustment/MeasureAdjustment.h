/*��������˺��������__declspec(dllexport)���η�����δ����������__declspec(dllimport)*/
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


//-----------ˮ׼��ƽ��-----------
//��λ���̣߳��ף����룺ǧ��
//���룺
//		Lnumber:�۲�߲�����
//		Pnumber���ܵ���
//		knPnumber����֪����
//		Sigma:��ǰ��λȨ�����,ÿǧ�׸߲���������Ϊ��λ��
// 
//		Pname:�������飬���е����洢�������У����������е�������Ӧ��ʵ�����
//		Height����֪�߳�ֵ���飬����������Ӧ�ĵ�����Ŷ�Ӧ,�ò���Ҳ���������
//		StartP:�۲�����  �߲����ţ�����洢������Ӧ����ţ�����
//		EndP:�۲�����  �߲��յ�ţ�����洢������Ӧ����ţ�����
//		L,P:�۲�����  �߲�ֵL���飬·�߳���P���飬������������StartP��EndP����������Ӧ
//����ֵ��
//      pvv:ƽ����pvv
//      mu:���λȨ�����
// 
//      Height���߳�ƽ��ֵ�����飩����ź͵������Ӧ
//      dx:�̸߳����������飩����ź͵������Ӧ��Height-dx��Ϊ���Ƹ߳�
//      m0:�߳��������飩����ź͵������Ӧ
// 
//      V:�߲�����������飩����ź͹۲�������Ŷ�Ӧ��L+V�Ǹ߲�ƽ��ֵ
//      P:�۲�߲��Ȩ�����飩
//      m1:�߲��������飩����ź͹۲�������Ŷ�Ӧ
MEASUREADJUSTMENT_API int LevelAD(const int Lnumber, const int Pnumber, const int knPnumber, const double Sigma, char** Pname, double* Height, int* StartP, int* EndP, double* L, double* P, double& pvv, double& mu, double* dx, double* m0, double* V, double* m1);


//ˮ׼������·�߱պϲ����closure error
//��λ���̣߳��ף����룺ǧ��
//���룺
//		Lnumber:�۲�߲�����
//		Pnumber���ܵ���
//		knPnumber����֪����
//		Sigma:��ǰ��λȨ�����
// 
//		Pname:�������飬���е����洢�������У����������е�������Ӧ��ʵ�����
//		Height����֪�߳�ֵ���飬����������Ӧ�ĵ�����Ŷ�Ӧ,�ò���Ҳ���������
//		StartP:�۲�����  �߲����ţ�����洢������Ӧ����ţ�����
//		EndP:�۲�����  �߲��յ�ţ�����洢������Ӧ����ţ�����
//		L,P:�۲�����  �߲�ֵL���飬·�߳���P���飬������������StartP��EndP����������Ӧ
//����ֵ��
//      line_name:����·���������飬����line_name[1]ΪA B C D F�����������֮��Ϊ�ո�
//      line_L:����·�߳�������
//      line_w:��Ӧ����·�ߵıպϲ�����
//      line_limit:��Ӧ����·�ߵ��޲�����
//      loop_name:���պ�·���������飬��loop_name[1]Ϊ A B C A,���������֮��Ϊ�ո�
//      loop_L:�պϻ�·�߳���
//      loop_w:�պϲ�
//      loop_limit:�޲�
MEASUREADJUSTMENT_API int LevelCE(const int Lnumber, const int Pnumber, const int knPnumber, const double Sigma, char** Pname, double* Height, int* StartP, int* EndP, double* L, double* P,
	std::vector<std::string>& line_name, std::vector<double>& line_L, std::vector<double>& line_w, std::vector<double>& line_limit,
	std::vector<std::string>& loop_name, std::vector<double>& loop_L, std::vector<double>& loop_w, std::vector<double>& loop_limit);

//ˮƽ����С����ƽ��(���굥λΪ��)
// ��λ���ǶȺͷ���۲�����Ϊ���ȣ��߳��۲ⵥλΪ��
// ���룺
//      Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:�ܵ�������֪����������������������ֵ�������߳��۲�ֵ��������λ�ǹ۲�ֵ����
//      ma,mS1,mS2,mT:����ֵ�����(����Ϊ��λ)���߳��̶�������Ϊ��λ�����߳��ı����������٣�����λ����������Ϊ��λ��
//      Pname:�������飬��СΪ���е��С
//      dir0,dir1,dir2,L:�洢����ֵ�۲����ݣ�dir0���������飬��С����������+1���洢���Ϊÿ�鷽�����ݵ�һ�й۲�ֵ��L�д洢����ţ�
//                      dir1��վ��ã���СΪ������������dir2��׼��ţ���СΪ����ֵ������L�洢����ֵ����С����ֵ������
//      S_dir1,S_dir2,S_L���洢�߳��۲����ݣ���С��Ϊ�߳��۲�ֵ������S_dir1��վ��ţ�S_dir2��׼��ţ�S_L�߳��۲�ֵ��
//      T_dir1,T_dir2,T_L:�洢��λ�ǹ۲����ݣ���С��Ϊ��λ�ǹ۲�������T_dir1��վ��ţ�T_dir2��׼��ţ�T_L��λ�ǹ۲�ֵ��
//      b_manuAP,b_manu_SP,b_manu_TP:�Ƿ����뷽��ֵ���߳�����λ��Ȩ��
//      manuAP,manuSP,manuTP:����ķ���ֵ���߳�����λ��Ȩ��
// //   LS_max:��С����ƽ����ѭ���в��޲��λ���ף�Ĭ��0.01
// 
// ���أ�V,S_V,T_V�ֱ��Ƿ���ֵ���߳����ͷ�λ�����
//      S_SVsΪS+V�Ǳ߳�ƽ��ֵ��T_TV��T+V�Ƿ�λ�ǹ۲�ƽ��ֵ
//      XY��mx,my,M�ǵ�����ƽ��ֵ�������
//      mu����ƽ���
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

//ˮƽ��Helmert����������ƣ���С����ƽ�
// ��λ���ǶȺͷ���۲�����Ϊ���ȣ��߳��۲ⵥλΪ��
// ���룺
//      Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:�ܵ�������֪����������������������ֵ�������߳��۲�ֵ��������λ�ǹ۲�ֵ����
//      ma,mS1,mS2,mT:����ֵ�����(����Ϊ��λ)���߳��̶�������Ϊ��λ�����߳��ı����������٣�����λ����������Ϊ��λ��
//      Pname:�������飬��СΪ���е��С
//      dir0,dir1,dir2,L:�洢����ֵ�۲����ݣ�dir0���������飬��С����������+1���洢���Ϊÿ�鷽�����ݵ�һ�й۲�ֵ��L�д洢����ţ�
//                      dir1��վ��ã���СΪ������������dir2��׼��ţ���СΪ����ֵ������L�洢����ֵ����С����ֵ������
//      S_dir1,S_dir2,S_L���洢�߳��۲����ݣ���С��Ϊ�߳��۲�ֵ������S_dir1��վ��ţ�S_dir2��׼��ţ�S_L�߳��۲�ֵ��
//      T_dir1,T_dir2,T_L:�洢��λ�ǹ۲����ݣ���С��Ϊ��λ�ǹ۲�������T_dir1��վ��ţ�T_dir2��׼��ţ�T_L��λ�ǹ۲�ֵ��
// 
// ���أ�V,S_V,T_V�ֱ��Ƿ���ֵ���߳����ͷ�λ�����
//      S_SVsΪS+V�Ǳ߳�ƽ��ֵ��T_TV��T+V�Ƿ�λ�ǹ۲�ƽ��ֵ
//      XY��mx,my,M�ǵ�����ƽ��ֵ�������
//      mu����ƽ���
//      HelCircle:Helmertѭ��������,Ĭ��10��
//
MEASUREADJUSTMENT_API int PlaneNetHelmertAd(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L, double* V, double* R_T, double* R_mT, double* R_S, double* R_ms,
    int* S_dir1, int* S_dir2, double* S_L, double* S_V, double* S_T, double* S_mT, double* S_SVs, double* S_ms,
    int* T_dir1, int* T_dir2, double* T_L, double* T_V, double* T_TV, double* T_mT, double* T_S, double* T_ms,
    double* XY, double* mx, double* my, double* M,
    double& mu, double& m1, double& m2,int HelCircle=10);

//ˮƽ��������ƽ��
// ��λ���ǶȺͷ���۲�����Ϊ���ȣ��߳��۲ⵥλΪ��
// ���룺
//      Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:�ܵ�������֪����������������������ֵ�������߳��۲�ֵ��������λ�ǹ۲�ֵ����
//      ma,mS1,mS2,mT:����ֵ�����(����Ϊ��λ)���߳��̶�������Ϊ��λ�����߳��ı����������٣�����λ����������Ϊ��λ��
//      Pname:�������飬��СΪ���е��С
//      dir0,dir1,dir2,L:�洢����ֵ�۲����ݣ�dir0���������飬��С����������+1���洢���Ϊÿ�鷽�����ݵ�һ�й۲�ֵ��L�д洢����ţ�
//                      dir1��վ��ã���СΪ������������dir2��׼��ţ���СΪ����ֵ������L�洢����ֵ����С����ֵ������
//      S_dir1,S_dir2,S_L���洢�߳��۲����ݣ���С��Ϊ�߳��۲�ֵ������S_dir1��վ��ţ�S_dir2��׼��ţ�S_L�߳��۲�ֵ��
//      T_dir1,T_dir2,T_L:�洢��λ�ǹ۲����ݣ���С��Ϊ��λ�ǹ۲�������T_dir1��վ��ţ�T_dir2��׼��ţ�T_L��λ�ǹ۲�ֵ��
// 
// ���أ�V,S_V,T_V�ֱ��Ƿ���ֵ���߳����ͷ�λ�����
//      S_SVsΪS+V�Ǳ߳�ƽ��ֵ��T_TV��T+V�Ƿ�λ�ǹ۲�ƽ��ֵ
//      XY��mx,my,M�ǵ�����ƽ��ֵ�������
//      mu����ƽ���
//
MEASUREADJUSTMENT_API int PlaneNetFree(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L, double* V, double* R_T, double* R_mT, double* R_S, double* R_ms,
    int* S_dir1, int* S_dir2, double* S_L, double* S_V, double* S_T, double* S_mT, double* S_SVs, double* S_ms,
    int* T_dir1, int* T_dir2, double* T_L, double* T_V, double* T_TV, double* T_mT, double* T_S, double* T_ms,
    double* XY, double* mx, double* my, double* M,
    double& mu);


//ˮƽ����ֵ����
// ��λ���ǶȺͷ���۲�����Ϊ���ȣ��߳��۲ⵥλΪ��
// ���룺
//      Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:�ܵ�������֪����������������������ֵ�������߳��۲�ֵ��������λ�ǹ۲�ֵ����
//      ma,mS1,mS2,mT:����ֵ�����(����Ϊ��λ)���߳��̶�������Ϊ��λ�����߳��ı����������٣�����λ����������Ϊ��λ��
//      Pname:�������飬��СΪ���е��С
//      dir0,dir1,dir2,L:�洢����ֵ�۲����ݣ�dir0���������飬��С����������+1���洢���Ϊÿ�鷽�����ݵ�һ�й۲�ֵ��L�д洢����ţ�
//                      dir1��վ��ã���СΪ������������dir2��׼��ţ���СΪ����ֵ������L�洢����ֵ����С����ֵ������
//      S_dir1,S_dir2,S_L���洢�߳��۲����ݣ���С��Ϊ�߳��۲�ֵ������S_dir1��վ��ţ�S_dir2��׼��ţ�S_L�߳��۲�ֵ��
//      T_dir1,T_dir2,T_L:�洢��λ�ǹ۲����ݣ���С��Ϊ��λ�ǹ۲�������T_dir1��վ��ţ�T_dir2��׼��ţ�T_L��λ�ǹ۲�ֵ��
// 
// ���أ�
//      XY:����ĳ�ʼ������
//
MEASUREADJUSTMENT_API int PlaneNetinitPoint(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L,
    int* S_dir1, int* S_dir2, double* S_L,
    int* T_dir1, int* T_dir2, double* T_L,
    double* XY);



//ˮƽ����������ֵ����
// ��λ���ǶȺͷ���۲�����Ϊ���ȣ��߳��۲ⵥλΪ��
// ���룺
//      Pnumber,knPnumber,Lnumber,Nnumber,Snumber,Tnumber:�ܵ�������֪����������������������ֵ�������߳��۲�ֵ��������λ�ǹ۲�ֵ����
//      ma,mS1,mS2,mT:����ֵ�����(����Ϊ��λ)���߳��̶�������Ϊ��λ�����߳��ı����������٣�����λ����������Ϊ��λ��
//      Pname:�������飬��СΪ���е��С
//      dir0,dir1,dir2,L:�洢����ֵ�۲����ݣ�dir0���������飬��С����������+1���洢���Ϊÿ�鷽�����ݵ�һ�й۲�ֵ��L�д洢����ţ�
//                      dir1��վ��ã���СΪ������������dir2��׼��ţ���СΪ����ֵ������L�洢����ֵ����С����ֵ������
//      S_dir1,S_dir2,S_L���洢�߳��۲����ݣ���С��Ϊ�߳��۲�ֵ������S_dir1��վ��ţ�S_dir2��׼��ţ�S_L�߳��۲�ֵ��
//      T_dir1,T_dir2,T_L:�洢��λ�ǹ۲����ݣ���С��Ϊ��λ�ǹ۲�������T_dir1��վ��ţ�T_dir2��׼��ţ�T_L��λ�ǹ۲�ֵ��
// 
// ���أ�
//      XY:����ĳ�ʼ������
//
MEASUREADJUSTMENT_API int PlaneNetFreeinitPoint(const int Pnumber, const int knPnumber, const int Lnumber, const int Nnumber, const int Snumber, const int Tnumber, const double ma, const double mS1, const double mS2, const double mT, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L,
    int* S_dir1, int* S_dir2, double* S_L,
    int* T_dir1, int* T_dir2, double* T_L,
    double* XY);


//ˮƽ�����պϲ����
//����ͬ����ͬ���֣�����ָ��һ��pair�飬������ɻ��ĵ�ͻ��պϲ�պϲλ�ǻ���
MEASUREADJUSTMENT_API std::vector<std::pair<std::vector<std::string>, float>> PlaneCE(const int Lnumber, const int Nnumber, char** Pname,
    int* dir0, int* dir1, int* dir2, double* L);