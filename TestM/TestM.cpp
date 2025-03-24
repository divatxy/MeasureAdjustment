// TestM.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include "../MeasureAdjustment/MeasureAdjustment.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>

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

void testlevelad()
{
    //FILE* fp = fopen("算例\\最小二乘平差\\data.txt", "r");
    std::ifstream file("..\\算例\\最小二乘平差\\data.txt");
    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << std::endl;
        return; // 如果文件打开失败，退出程序
    }
    int Lnumber(0), Pnumber(0), knPnumber(0);
    double sigma(0);

    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> Lnumber >> Pnumber >> knPnumber >> sigma;
    //fscanf(fp, "%d%d%d%", Lnumber, Pnumber, knPnumber);
    double* height = new double[Pnumber];
    int* startp = new int[Lnumber];
    int* endp = new int[Lnumber];
    double* L = new double[Lnumber];
    double* p = new double[Lnumber];
    char** pname = new char* [Pnumber];
    for (size_t i = 0; i < Pnumber; i++)
    {
        pname[i] = NULL;
    }
    int inum(0);
    while (std::getline(file, line))
    {
        if (line.empty())// std::string的empty()方法检查字符串是否为空
        {
            continue;
        }
        else
        {
            std::istringstream iss1(line);
            //读取已知高程
            if (inum < knPnumber)
            {
                char* buffer = new char[100];
                iss1 >> buffer;
                int index;
                for (size_t i = 0; i < Pnumber; i++)
                {
                    if (pname[i] != NULL)
                    {
                        //将待查点名与已经存入点名数组的点名比较
                        if (strcmp(buffer, pname[i]) == 0)
                        {
                            index = i;
                            break;
                        }
                    }
                    else
                    {
                        //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                        int len = strlen(buffer);
                        pname[i] = new char[len + 1];
                        strcpy_s(pname[i], len + 1, buffer);
                        index = i;
                        break;
                    }
                }
                double d1;
                iss1 >> d1;
                height[index] = d1;
            }
            if ((inum + 1) == knPnumber)
            {
                break;
            }
            inum++;
        }
    }
    //读取观测数据
    inum = 0;
    while (std::getline(file, line))
    {
        if (!line.empty())// std::string的empty()方法检查字符串是否为空
        {
            std::istringstream iss1(line);
            if (inum < Lnumber)
            {
                char* buffer = new char[100];
                iss1 >> buffer;//读取高程起点
                int index;
                for (size_t i = 0; i < Pnumber; i++)
                {
                    if (pname[i] != NULL)
                    {
                        //将待查点名与已经存入点名数组的点名比较
                        if (strcmp(buffer, pname[i]) == 0) { index = i; break; }
                    }
                    else
                    {
                        //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                        int len = strlen(buffer);
                        pname[i] = new char[len + 1];
                        strcpy_s(pname[i], len + 1, buffer);
                        index = i;
                        break;
                    }
                }
                startp[inum] = index;

                iss1 >> buffer;//读取高程终点
                //int index;
                for (size_t i = 0; i < Pnumber; i++)
                {
                    if (pname[i] != NULL)
                    {
                        //将待查点名与已经存入点名数组的点名比较
                        if (strcmp(buffer, pname[i]) == 0) { index = i; break; }
                    }
                    else
                    {
                        //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                        int len = strlen(buffer);
                        pname[i] = new char[len + 1];
                        strcpy_s(pname[i], len + 1, buffer);
                        index = i;
                        break;
                    }
                }
                endp[inum] = index;
                //读取高差值与路线长度
                double d1, d2;
                iss1 >> d1 >> d2;
                L[inum] = d1;
                p[inum] = d2;
            }
            if ((inum + 1) == Lnumber)
            {
                break;
            }
            inum++;
        }
    }
    //数据读取完毕进行处理
    double pvv, mu;
    double* dx = new double[Pnumber];
    double* m0 = new double[Pnumber];
    double* v = new double[Lnumber];
    double* m1 = new double[Lnumber];
    LevelAD(Lnumber, Pnumber, knPnumber, sigma, pname, height, startp, endp, L, p, pvv, mu, dx, m0, v, m1);
    std::cout << "pvv:" << pvv << "\t mu:" << mu << std::endl;
    std::cout << "==== 高程平差值及其精度 ====" << std::endl;
    std::cout << "点名   近似高程   改正数   高程平差值  中误差" << std::endl;
    for (size_t i = 0; i < Pnumber; i++)
    {
        std::cout << pname[i] << "\t" << height[i] - dx[i] << "\t" << dx[i] << "\t" << height[i] << "\t" << m0[i] << std::endl;
    }
    std::cout << "==== 观测值平差值及其精度 ====" << std::endl;
    std::cout << "No. 起点  终点  观测高差    ｖ   高差平差值   观测权   中误差" << std::endl;
    for (size_t i = 0; i < Lnumber; i++)
    {
        int k1 = startp[i];
        int k2 = endp[i];
        std::cout << i + 1 << "\t" << pname[k1] << "\t" << pname[k2] << "\t" << L[i] << "\t" << v[i] << "\t" << L[i] + v[i] << "\t" << p[i] << "\t" << m1[i] << std::endl;
    }
}

void testlevelce()
{
    std::ifstream file("..\\算例\\闭合差计算\\Data.txt");
    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << std::endl;
        return; // 如果文件打开失败，退出程序
    }
    int Lnumber(0), Pnumber(0), knPnumber(0);
    double sigma(0);

    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> Lnumber >> Pnumber >> knPnumber >> sigma;
    //fscanf(fp, "%d%d%d%", Lnumber, Pnumber, knPnumber);
    double* height = new double[Pnumber];
    int* startp = new int[Lnumber];
    int* endp = new int[Lnumber];
    double* L = new double[Lnumber];
    double* p = new double[Lnumber];
    char** pname = new char* [Pnumber];
    for (size_t i = 0; i < Pnumber; i++)
    {
        pname[i] = NULL;
    }
    int inum(0);
    while (std::getline(file, line))
    {
        if (line.empty())// std::string的empty()方法检查字符串是否为空
        {
            continue;
        }
        else
        {
            std::istringstream iss1(line);
            //读取已知高程
            if (inum < knPnumber)
            {
                char* buffer = new char[100];
                iss1 >> buffer;
                int index;
                for (size_t i = 0; i < Pnumber; i++)
                {
                    if (pname[i] != NULL)
                    {
                        //将待查点名与已经存入点名数组的点名比较
                        if (strcmp(buffer, pname[i]) == 0)
                        {
                            index = i;
                            break;
                        }
                    }
                    else
                    {
                        //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                        int len = strlen(buffer);
                        pname[i] = new char[len + 1];
                        strcpy_s(pname[i], len + 1, buffer);
                        index = i;
                        break;
                    }
                }
                double d1;
                iss1 >> d1;
                height[index] = d1;
            }
            if ((inum + 1) == knPnumber)
            {
                break;
            }
            inum++;
        }
    }
    //读取观测数据
    inum = 0;
    while (std::getline(file, line))
    {
        if (!line.empty())// std::string的empty()方法检查字符串是否为空
        {
            std::istringstream iss1(line);
            if (inum < Lnumber)
            {
                char* buffer = new char[100];
                iss1 >> buffer;//读取高程起点
                int index;
                for (size_t i = 0; i < Pnumber; i++)
                {
                    if (pname[i] != NULL)
                    {
                        //将待查点名与已经存入点名数组的点名比较
                        if (strcmp(buffer, pname[i]) == 0) { index = i; break; }
                    }
                    else
                    {
                        //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                        int len = strlen(buffer);
                        pname[i] = new char[len + 1];
                        strcpy_s(pname[i], len + 1, buffer);
                        index = i;
                        break;
                    }
                }
                startp[inum] = index;

                iss1 >> buffer;//读取高程终点
                //int index;
                for (size_t i = 0; i < Pnumber; i++)
                {
                    if (pname[i] != NULL)
                    {
                        //将待查点名与已经存入点名数组的点名比较
                        if (strcmp(buffer, pname[i]) == 0) { index = i; break; }
                    }
                    else
                    {
                        //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                        int len = strlen(buffer);
                        pname[i] = new char[len + 1];
                        strcpy_s(pname[i], len + 1, buffer);
                        index = i;
                        break;
                    }
                }
                endp[inum] = index;
                //读取高差值与路线长度
                double d1, d2;
                iss1 >> d1 >> d2;
                L[inum] = d1;
                p[inum] = d2;
            }
            if ((inum + 1) == Lnumber)
            {
                break;
            }
            inum++;
        }
    }
    //数据读取完毕进行处理
    std::vector<std::string> li_name, lo_name;
    std::vector<double> li_L, lo_L;
    std::vector<double> li_w, lo_w;
    std::vector<double> li_lim, lo_lim;
    LevelCE(Lnumber, Pnumber, knPnumber, sigma, pname, height, startp, endp, L, p, li_name, li_L, li_w, li_lim, lo_name, lo_L, lo_w, lo_lim);
    int li_num = li_name.size();
    int lo_num = lo_name.size();
    std::cout << "\n==== 附合路线闭合差 ====" << std::endl;
    std::cout << "序号   路线   路线长度    闭合差  限差" << std::endl;
    for (size_t i = 0; i < li_num; i++)
    {
        std::cout << i + 1 << "\t" << li_name[i] << "\t" << li_L[i] << "\t" << li_w[i] << "\t" << li_lim[i] << std::endl;
    }
    std::cout << "\n==== 闭合环路线闭合差 ====" << std::endl;
    std::cout << "序号   路线   路线长度    闭合差  限差" << std::endl;
    for (size_t i = 0; i < lo_num; i++)
    {
        std::cout << i + 1 << "\t" << lo_name[i] << "\t" << lo_L[i] << "\t" << lo_w[i] << "\t" << lo_lim[i] << std::endl;
    }
}

void testPlaneAd() 
{
    std::ifstream file("..\\算例\\最小二乘平差\\Data1.txt");
    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << std::endl;
        return; // 如果文件打开失败，退出程序
    }
    int Pnumber, knPnumber, Lnumber, Nnumber, Snumber, Tnumber;
    double ma, ms1, ms2, mt;
    
    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> Pnumber >> knPnumber >> Lnumber >> Nnumber >> Snumber >> Tnumber;

    std::getline(file, line);//这一行是空行
    std::getline(file, line);
    iss.str(line);
    iss >> ma >> ms1 >> ms2 >> mt;
    
    char** pname = new char* [Pnumber];//点名
    for (size_t i = 0; i < Pnumber; i++)
    {
        pname[i] = NULL;
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
    if (Lnumber>0)
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
    int inum(0);
    while (std::getline(file, line))//读取已知点坐标
    {
        if (line.empty())// std::string的empty()方法检查字符串是否为空
        {
            continue;
        }
        else
        {
            std::istringstream iss1(line);
            if (inum < knPnumber)
            {
                char* buffer = new char[100];
                iss1 >> buffer;
                int index;
                for (size_t i = 0; i < Pnumber; i++)
                {
                    if (pname[i] != NULL)
                    {
                        //将待查点名与已经存入点名数组的点名比较
                        if (strcmp(buffer, pname[i]) == 0)
                        {
                            index = i;
                            break;
                        }
                    }
                    else
                    {
                        //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                        int len = strlen(buffer);
                        pname[i] = new char[len + 1];
                        strcpy_s(pname[i], len + 1, buffer);
                        index = i;
                        break;
                    }
                }
                double x,y;
                iss1 >> x >> y;
                XY[2*index] = x;
                XY[2 * index + 1] = y;
            }
            if ((inum + 1) == knPnumber)
            {
                break;
            }
            inum++;
        }
    }
    inum = 0;
    
    dir0[0] = 0;
    for (size_t i = 0; i <= Lnumber-1; i++)//读取方向值
    {
        std::getline(file, line);//空行
        std::getline(file, line);
        std::istringstream iss1(line);
        char* buffer = new char[100];
        int ni = 0;
        iss1 >> buffer >> ni;
        //存储测站点名
        int index=-1;
        for (size_t j = 0; j < Pnumber; j++)
        {
            if (pname[j] != NULL)
            {
                //将待查点名与已经存入点名数组的点名比较
                if (strcmp(buffer, pname[j]) == 0)
                {
                    index = j;
                    break;
                }
            }
            else
            {
                //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                int len = strlen(buffer);
                pname[j] = new char[len + 1];
                strcpy_s(pname[j], len + 1, buffer);
                index = j;
                break;
            }
			if (j == (Pnumber-1) && index == -1)
            {
                std::cout << "error" << std::endl;
                return ;
            }
        }
        dir1[i] = index;
        dir0[i + 1] = dir0[i] + ni;
        for (size_t ii = dir0[i]; ii < dir0[i+1]; ii++)
        {
            std::getline(file, line);
            iss1.clear();
            iss1.str(line);
            iss1 >> buffer;
            for (size_t j = 0; j < Pnumber; j++)
            {
                if (pname[j] != NULL)
                {
                    //将待查点名与已经存入点名数组的点名比较
                    if (strcmp(buffer, pname[j]) == 0)
                    {
                        index = j;
                        break;
                    }
                }
                else
                {
                    //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                    int len = strlen(buffer);
                    pname[j] = new char[len + 1];
                    strcpy_s(pname[j], len + 1, buffer);
                    index = j;
                    break;
                }
            }
            dir2[ii] = index;//照准点号
            double Lvalue;
            iss1 >> Lvalue;
            L[ii] = dms_rad(Lvalue);
        }
    }
    //读取边长
    std::getline(file, line);//空行
    for (size_t i = 0; i < Snumber; i++)
    {
        std::getline(file, line);
        std::istringstream iss1(line);
        char* buffer = new char[100];
        iss1 >> buffer;
        //存储测站点名
        int index;
        for (size_t j = 0; j < Pnumber; j++)
        {
            if (pname[j] != NULL)
            {
                //将待查点名与已经存入点名数组的点名比较
                if (strcmp(buffer, pname[j]) == 0)
                {
                    index = j;
                    break;
                }
            }
            else
            {
                //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                int len = strlen(buffer);
                pname[j] = new char[len + 1];
                strcpy_s(pname[j], len + 1, buffer);
                index = j;
                break;
            }
        }
        S_dir1[i] = index;
        iss1 >> buffer;
        //存储测站点名
        for (size_t j = 0; j < Pnumber; j++)
        {
            if (pname[j] != NULL)
            {
                //将待查点名与已经存入点名数组的点名比较
                if (strcmp(buffer, pname[j]) == 0)
                {
                    index = j;
                    break;
                }
            }
            else
            {
                //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                int len = strlen(buffer);
                pname[j] = new char[len + 1];
                strcpy_s(pname[j], len + 1, buffer);
                index = j;
                break;
            }
        }
        S_dir2[i] = index;
        iss1 >> S_L[i];
    }
    //读取方位角
    std::getline(file, line);//空行
    for (size_t i = 0; i < Tnumber; i++)
    {
        std::getline(file, line);
        std::istringstream iss1(line);
        char* buffer = new char[100];
        iss1 >> buffer;
        //存储测站点名
        int index;
        for (size_t j = 0; j < Pnumber; j++)
        {
            if (pname[j] != NULL)
            {
                //将待查点名与已经存入点名数组的点名比较
                if (strcmp(buffer, pname[j]) == 0)
                {
                    index = j;
                    break;
                }
            }
            else
            {
                //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                int len = strlen(buffer);
                pname[j] = new char[len + 1];
                strcpy_s(pname[j], len + 1, buffer);
                index = j;
                break;
            }
        }
        T_dir1[i] = index;
        iss1 >> buffer;
        //存储测站点名
        for (size_t j = 0; j < Pnumber; j++)
        {
            if (pname[j] != NULL)
            {
                //将待查点名与已经存入点名数组的点名比较
                if (strcmp(buffer, pname[j]) == 0)
                {
                    index = j;
                    break;
                }
            }
            else
            {
                //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                int len = strlen(buffer);
                pname[j] = new char[len + 1];
                strcpy_s(pname[j], len + 1, buffer);
                index = j;
                break;
            }
        }
        T_dir2[i] = index;
        iss1 >> T_L[i];
        T_L[i] = dms_rad(T_L[i]);
    }
    double mu = 0;//平差总精度
    PlaneNetAD(Pnumber, knPnumber, Lnumber, Nnumber, Snumber, Tnumber, ma, ms1, ms2, mt, pname,
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
        for (size_t j = dir0[i]; j < dir0[i+1]; j++)
        {
            if (j==dir0[i])
            {
                std::cout << "\n" << pname[dir1[i]];

            }
            else
            {
                //std::cout << "\n";
            }
            std::cout << "\t" << pname[dir2[j]] << "\t" << rad_dms(L[j]) << "\t" << V[j] << "\t" << rad_dms(R_T[j]) << "\t" << R_mt[j] << "\t" << R_S[j] << "\t" << R_ms[j] << std::endl;
        }
    }

    std::cout << "\n" << std::endl;
    std::cout << "------边长观测值平差成果-----" << std::endl;
    std::cout << "P1 \t P2 \t S \t Vs \t T \t mT \t S+Vs \t ms" << std::endl;
    for (size_t i = 0; i < Snumber; i++)
    {
        std::cout << pname[S_dir1[i]] << "\t" << pname[S_dir2[i]] << "\t" << S_L[i] << "\t" << S_V[i] << "\t" << rad_dms(S_T[i]) << "\t" << S_mt[i] << "\t" << S_S[i] << "\t" << S_ms[i] << std::endl;
    }

    std::cout << "\n" << std::endl;
    std::cout << "------方位角平差成果-----" << std::endl;
    std::cout << "P1 \t P2 \t S \t Vs \t T \t mT \t S+Vs \t ms" << std::endl;
    for (size_t i = 0; i < Tnumber; i++)
    {
        std::cout << pname[T_dir1[i]] << "\t" << pname[T_dir2[i]] << "\t" << T_L[i] << "\t" << T_V[i] << "\t" << rad_dms(T_T[i]) << "\t" << T_mt[i] << "\t" << T_S[i] << "\t" << T_ms[i] << std::endl;
    }

}

//导线网自由网平差
void testPlaneAdFree()
{
    std::ifstream file("..\\算例\\自由网平差\\Data.txt");
    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << std::endl;
        return; // 如果文件打开失败，退出程序
    }
    int Pnumber, knPnumber, Lnumber, Nnumber, Snumber, Tnumber;
    double ma, ms1, ms2, mt;

    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> Pnumber >> knPnumber >> Lnumber >> Nnumber >> Snumber >> Tnumber;

    std::getline(file, line);//这一行是空行
    std::getline(file, line);
    iss.str(line);
    iss >> ma >> ms1 >> ms2 >> mt;

    char** pname = new char* [Pnumber];//点名
    for (size_t i = 0; i < Pnumber; i++)
    {
        pname[i] = NULL;
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
    int inum(0);

    dir0[0] = 0;
    for (size_t i = 0; i <= Lnumber - 1; i++)//读取方向值
    {
        std::getline(file, line);//空行
        std::getline(file, line);
        std::istringstream iss1(line);
        char* buffer = new char[100];
        int ni = 0;
        iss1 >> buffer >> ni;
        //存储测站点名
        int index = -1;
        for (size_t j = 0; j < Pnumber; j++)
        {
            if (pname[j] != NULL)
            {
                //将待查点名与已经存入点名数组的点名比较
                if (strcmp(buffer, pname[j]) == 0)
                {
                    index = j;
                    break;
                }
            }
            else
            {
                //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                int len = strlen(buffer);
                pname[j] = new char[len + 1];
                strcpy_s(pname[j], len + 1, buffer);
                index = j;
                break;
            }
            if (j == (Pnumber - 1) && index == -1)
            {
                std::cout << "error" << std::endl;
                return;
            }
        }
        dir1[i] = index;
        dir0[i + 1] = dir0[i] + ni;
        for (size_t ii = dir0[i]; ii < dir0[i + 1]; ii++)
        {
            std::getline(file, line);
            iss1.clear();
            iss1.str(line);
            iss1 >> buffer;
            for (size_t j = 0; j < Pnumber; j++)
            {
                if (pname[j] != NULL)
                {
                    //将待查点名与已经存入点名数组的点名比较
                    if (strcmp(buffer, pname[j]) == 0)
                    {
                        index = j;
                        break;
                    }
                }
                else
                {
                    //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                    int len = strlen(buffer);
                    pname[j] = new char[len + 1];
                    strcpy_s(pname[j], len + 1, buffer);
                    index = j;
                    break;
                }
            }
            dir2[ii] = index;//照准点号
            double Lvalue;
            iss1 >> Lvalue;
            L[ii] = dms_rad(Lvalue);
        }
    }
    //读取边长
    std::getline(file, line);//空行
    for (size_t i = 0; i < Snumber; i++)
    {
        std::getline(file, line);
        std::istringstream iss1(line);
        char* buffer = new char[100];
        iss1 >> buffer;
        //存储测站点名
        int index;
        for (size_t j = 0; j < Pnumber; j++)
        {
            if (pname[j] != NULL)
            {
                //将待查点名与已经存入点名数组的点名比较
                if (strcmp(buffer, pname[j]) == 0)
                {
                    index = j;
                    break;
                }
            }
            else
            {
                //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                int len = strlen(buffer);
                pname[j] = new char[len + 1];
                strcpy_s(pname[j], len + 1, buffer);
                index = j;
                break;
            }
        }
        S_dir1[i] = index;
        iss1 >> buffer;
        //存储测站点名
        for (size_t j = 0; j < Pnumber; j++)
        {
            if (pname[j] != NULL)
            {
                //将待查点名与已经存入点名数组的点名比较
                if (strcmp(buffer, pname[j]) == 0)
                {
                    index = j;
                    break;
                }
            }
            else
            {
                //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                int len = strlen(buffer);
                pname[j] = new char[len + 1];
                strcpy_s(pname[j], len + 1, buffer);
                index = j;
                break;
            }
        }
        S_dir2[i] = index;
        iss1 >> S_L[i];
    }
    //读取方位角
    std::getline(file, line);//空行
    for (size_t i = 0; i < Tnumber; i++)
    {
        std::getline(file, line);
        std::istringstream iss1(line);
        char* buffer = new char[100];
        iss1 >> buffer;
        //存储测站点名
        int index;
        for (size_t j = 0; j < Pnumber; j++)
        {
            if (pname[j] != NULL)
            {
                //将待查点名与已经存入点名数组的点名比较
                if (strcmp(buffer, pname[j]) == 0)
                {
                    index = j;
                    break;
                }
            }
            else
            {
                //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                int len = strlen(buffer);
                pname[j] = new char[len + 1];
                strcpy_s(pname[j], len + 1, buffer);
                index = j;
                break;
            }
        }
        T_dir1[i] = index;
        iss1 >> buffer;
        //存储测站点名
        for (size_t j = 0; j < Pnumber; j++)
        {
            if (pname[j] != NULL)
            {
                //将待查点名与已经存入点名数组的点名比较
                if (strcmp(buffer, pname[j]) == 0)
                {
                    index = j;
                    break;
                }
            }
            else
            {
                //待查点名是一个新的点名，将新点名的地址放到Pname数组中
                int len = strlen(buffer);
                pname[j] = new char[len + 1];
                strcpy_s(pname[j], len + 1, buffer);
                index = j;
                break;
            }
        }
        T_dir2[i] = index;
        iss1 >> T_L[i];
        T_L[i] = dms_rad(T_L[i]);
    }
    double mu = 0;//平差总精度
    PlaneNetFree(Pnumber, knPnumber, Lnumber, Nnumber, Snumber, Tnumber, ma, ms1, ms2, mt, pname,
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
            std::cout << "\t" << pname[dir2[j]] << "\t" << rad_dms(L[j]) << "\t" << V[j] << "\t" << rad_dms(R_T[j]) << "\t" << R_mt[j] << "\t" << R_S[j] << "\t" << R_ms[j] << std::endl;
        }
    }

    std::cout << "\n" << std::endl;
    std::cout << "------边长观测值平差成果-----" << std::endl;
    std::cout << "P1 \t P2 \t S \t Vs \t T \t mT \t S+Vs \t ms" << std::endl;
    for (size_t i = 0; i < Snumber; i++)
    {
        std::cout << pname[S_dir1[i]] << "\t" << pname[S_dir2[i]] << "\t" << S_L[i] << "\t" << S_V[i] << "\t" << rad_dms(S_T[i]) << "\t" << S_mt[i] << "\t" << S_S[i] << "\t" << S_ms[i] << std::endl;
    }

    std::cout << "\n" << std::endl;
    std::cout << "------方位角平差成果-----" << std::endl;
    std::cout << "P1 \t P2 \t S \t Vs \t T \t mT \t S+Vs \t ms" << std::endl;
    for (size_t i = 0; i < Tnumber; i++)
    {
        std::cout << pname[T_dir1[i]] << "\t" << pname[T_dir2[i]] << "\t" << T_L[i] << "\t" << T_V[i] << "\t" << rad_dms(T_T[i]) << "\t" << T_mt[i] << "\t" << T_S[i] << "\t" << T_ms[i] << std::endl;
    }

}

int main()
{
    //testlevelad();
    //testlevelce();
    testPlaneAd();
    //testPlaneAdFree();
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
