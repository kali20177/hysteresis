#include "pch.h"
#include <Windows.h>
#include <matplotlibcpp.h>
#include <fstream>

using namespace matplotlibcpp;

#define M_PI 3.14159265358979323846

double interval = 0.00002; // 50 hz 时间间隔 count = 1000

//double shape_control_fun(int data_num, int data_size, double* arr);

template <typename T>
std::vector<T>& operator+(std::vector<T>& v1, std::vector<T>& v2)
{
	v1.insert(v1.end(), v2.begin(), v2.end());
	return v1;
}


//产生两点之间序列
//double temp = 0.0;

int count = 1000; // 周期点数

std::vector<double> generate_vol_seq(const double& v1, const double& v2)
{
	std::vector<double> vol_seq(count, 0);

	for (int i = 0; i < count; ++i)
	{
		double vol;
		if (v2 >= v1)
			vol = fabs((v2 - v1) / 2) * sin(2 * M_PI * (i + 1) / (2 * count) - M_PI / 2) + fabs((v2 - v1) / 2) + min(v2, v1);
		else
			vol = fabs((v2 - v1) / 2) * sin(2 * M_PI * (i + 1) / (2 * count) + M_PI / 2) + fabs((v2 - v1) / 2) + min(v2, v1);

		vol_seq[i] = vol;
	}
	//temp = new_vol;
	return vol_seq;
}

std::vector<double> generate_disp_seq(const double& d1, const double& d2)
{
	std::vector<double> disp_seq(count, 0);

	for (int i = 0; i < count; ++i)
	{
		double disp;
		if (d2 >= d1)
			disp = fabs((d2 - d1) / 2) * sin(2 * M_PI * (i + 1) / (2 * count) - M_PI / 2) + fabs((d2 - d1) / 2) + min(d2, d1);
		else
			disp = fabs((d2 - d1) / 2) * sin(2 * M_PI * (i + 1) / (2 * count) + M_PI / 2) + fabs((d2 - d1) / 2) + min(d2, d1);

		disp_seq[i] = disp;
	}
	return disp_seq;
}

//double shape_control_fun(int data_num, int data_size, double* arr)
//{
//	double mid_pra = 0;
//	double part = (double)data_num / data_size;
//
//	if (part < 0.125)
//	{
//		mid_pra = *(arr + 3);
//	}
//	else if (part < 0.25)
//	{
//		mid_pra = *(arr + 4);
//	}
//	else if (part < 0.375)
//	{
//		mid_pra = *(arr + 5);
//	}
//	else if (part < 0.5)
//	{
//		mid_pra = *(arr + 6);
//	}
//	else if (part < 0.625)
//	{
//		mid_pra = *(arr + 7);
//	}
//	else if (part < 0.75)
//	{
//		mid_pra = *(arr + 8);
//	}
//	else if (part < 0.875)
//	{
//		mid_pra = *(arr + 9);
//	}
//	else
//	{
//		mid_pra = *(arr + 10);
//	}
//
//	return mid_pra;
//}

double shape_control_fun(int data_num, int data_size, const double* parameter)
{
	double mid_pra = 0;
	double part = (double)data_num / data_size;

	if (part < 0.125)
		mid_pra = *(parameter + 3);
	else if (part < 0.25)
		mid_pra = *(parameter + 4);
	else if (part < 0.375)
		mid_pra = *(parameter + 5);
	else if (part < 0.5)
		mid_pra = *(parameter + 6);
	else if (part < 0.625)
		mid_pra = *(parameter + 7);
	else if (part < 0.75)
		mid_pra = *(parameter + 8);
	else if (part < 0.875)
		mid_pra = *(parameter + 9);
	else
		mid_pra = *(parameter + 10);
	return mid_pra;
}

double shape_control_fun1(int data_num, int data_size, const double* parameter)
{
	double mid_pra = 0;
	double part = (double)data_num / data_size;
	if (part < 0.25)
		mid_pra = *(parameter + 3);
	else if (part < 0.5)
		mid_pra = *(parameter + 4);
	else if (part < 0.75)
		mid_pra = *(parameter + 5);
	else
		mid_pra = *(parameter + 6);
	return mid_pra;
}

double shape_control_fun2(int data_num, int data_size, const double* parameter)
{
	double mid_pra = 0;
	double part = (double)data_num / data_size;
	if (part < 0.25)
		mid_pra = *(parameter + 7);
	else if (part < 0.5)
		mid_pra = *(parameter + 8);
	else if (part < 0.75)
		mid_pra = *(parameter + 9);
	else
		mid_pra = *(parameter + 10);
	return mid_pra;
}

std::vector<double> feed_forward(const double* parameter, const std::vector<double>& input_displace)
{
	double k, a, y, eta = 0;
	//double try_interval =(double) interval;

	k = *parameter;
	a = *(parameter + 1);
	y = *(parameter + 2);
	eta = *(parameter + 11);

	double mid_prameter = 0;

	std::vector<double> u_ff(input_displace.size(), 0);  //前馈电压
	std::vector<double> d_uff(input_displace.size(), 0);  //前馈电压倒数
	std::vector<double> d_dis(input_displace.size(), 0);  //目标位移导数
	std::vector<double> mid_h(input_displace.size(), 0);  //中间变量

	for (int i = 0; i < input_displace.size() - 1; ++i)
	{
		d_dis[i] = ((input_displace[i + 1] - input_displace[i]) / interval);
	}

	u_ff[1] = 1 / k * (input_displace[1] + mid_h[0] + eta * d_dis[1]);

	for (int i = 2; i < input_displace.size(); ++i)
	{
		const int m = (i - 2) / count;
		const int n = (i - 2) % count;

		if (m % 2 == 0)
			mid_prameter = shape_control_fun1(n, count, parameter);
		else
			mid_prameter = shape_control_fun2(n, count, parameter);

		d_uff[i - 2] = (u_ff[i - 1] - u_ff[i - 2]) / interval;
		mid_h[i - 1] = mid_h[i - 2] + interval * d_uff[i - 2] * (a - fabs(mid_h[i - 2]) * (y + mid_prameter));
		u_ff[i] = 1 / k * (input_displace[i] + mid_h[i - 1] + eta * d_dis[i]);
	}
	return u_ff;
}

std::vector<double> predict_disp(const double* parameter, const std::vector<double>& input_vol)
{
	double k, a, y, eta = 0;
	double mid_prameter = 0;

	k = *parameter;
	a = *(parameter + 1);
	y = *(parameter + 2);
	eta = *(parameter + 11);

	std::vector<double> d_input(input_vol.size(), 0);
	std::vector<double> e_hystery(input_vol.size(), 0);
	std::vector<double> e_disp(input_vol.size(), 0);

	for (int i = 0; i < input_vol.size(); ++i)
	{
		if (i < input_vol.size() - 1)
			d_input[i] = ((input_vol[i + 1] - input_vol[i]) / interval);
		else
			d_input[i] = d_input[i - 1];
	}

	for (int i = 0; i < input_vol.size() - 1; ++i)
	{
		const int m = i / count;
		const int n = i % count;

		if (m % 2 == 0)
			mid_prameter = shape_control_fun1(n, count, parameter);
		else
			mid_prameter = shape_control_fun2(n, count, parameter);

		/*if (i < count + 1)
			mid_prameter = shape_control_fun(i, count, parameter);
		else if (i < count * 2)
			mid_prameter = shape_control_fun(i - count, count, parameter);
		else if (i < count * 3)
			mid_prameter = shape_control_fun(i - 2 * count, count, parameter);
		else if (i < count * 4)
			mid_prameter = shape_control_fun(i - 3 * count, count, parameter);
		else if (i < count * 5)
			mid_prameter = shape_control_fun(i - 4 * count, count, parameter);
		else if (i < count * 6)
			mid_prameter = shape_control_fun(i - 5 * count, count, parameter);
		else if (i < count * 7)
			mid_prameter = shape_control_fun(i - 6 * count, count, parameter);
		else if (i < count * 8)
			mid_prameter = shape_control_fun(i - 7 * count, count, parameter);
		else if (i < count * 9)
			mid_prameter = shape_control_fun(i - 8 * count, count, parameter);
		else if (i < count * 10)
			mid_prameter = shape_control_fun(i - 9 * count, count, parameter);*/

		e_hystery[i + 1] = e_hystery[i] + interval * d_input[i] * (a - fabs(e_hystery[i]) * (y + mid_prameter));
	}

	for (int i = 0; i < e_hystery.size() - 1; ++i)
	{
		e_disp[i + 1] = ((eta - interval) * e_disp[i] + (k * input_vol[i] - e_hystery[i]) * interval) / eta;
	}

	std::ofstream p;
	p.open("output.csv", std::ios::out | std::ios::trunc);                //打开文件路径

	p << "vol" << "," << "d_vol" << "," << "hystery" << "," << "disp" << "\n";

	for (int i = 0; i < input_vol.size(); i++)
	{
		p << input_vol[i] << "," << d_input[i] << "," << e_hystery[i] << "," << e_disp[i] << "\n";
	}
	p.close();

	return e_disp;
}


int main()
{
	/* 不同频率参数 */
	/*double parameter_0_5[12] = { 0.108237, 0.035945,0.044645 ,0.00938 , 0.025849  ,0.046276,-0.000929, -0.034877, -0.023104, 0.019958, 0.066842, 0.032767 };
	double parameter_1[12] = { 0.115663, 0.032714, 0.065055, -0.037197, -0.036097, -0.032535, -0.064296, -0.071229, -0.070129, 0.051244, 0.019675, 0.019869 };
	double parameter_10[12] = { 0.108122 ,0.033320	,0.002089 ,-0.015927 ,0.022808 ,0.038053 ,0.023244 ,-0.025724 ,0.002076 ,0.036451 ,0.057263 ,0.0012229 };
	double parameter_20[12] = { 0.107275 ,0.039021	,0.011291,-0.001252 ,0.024191 ,0.038566	,0.022611 ,-0.019637 ,0.011842 ,0.032542 ,0.053771 ,0.000715 };
	double parameter_40[12] = { 0.106873,0.032097,0.042661,-0.059413 ,-0.020232	,-0.002878,-0.020551,-0.065032 ,-0.043613,-0.007908,0.006631,0.00043 };


	double parameter_75[12] = { 0.103219,	0.031346,0.0236024,-0.0251676,-0.0045288,0.0192964,0.0026112,-0.032903,-0.016674,0.0208954,0.0223352,0.000307 };
	double parameter_80[12] = { 0.101919,	0.030607,	0.031162, -0.02816,	-0.002262,	0.021447,	-0.004282, -0.048626,	-0.008386,	0.005682,	0.021782,	0.000327 };
	double parameter_100[12] = { 0.103428,	0.036161,	0.021968, -0.016372,	0.012578,	0.022066,	0.012497, -0.024407,	0.001386,	0.02069,	0.039607,	0.000245 };*/

	double parameter_50[12] = { 0.106678, 0.033195, 0.040417, -0.048901, -0.013362,	0.003823, -0.012365, -0.062161, -0.030806, -0.006678, 0.005098, 0.000386 };

	/*std::vector<double> temp_vol1 = generate_vol_seq(0.0, 150.0);
	std::vector<double> temp_vol2 = generate_vol_seq(150.0, 0.0);
	std::vector<double> temp_vol3 = generate_vol_seq(0.0, 150.0);
	std::vector<double> temp_vol4 = generate_vol_seq(150.0, 20.0);
	std::vector<double> temp_vol5 = generate_vol_seq(20.0, 100.0);
	std::vector<double> temp_vol6 = generate_vol_seq(100.0, 70.0);
	std::vector<double> temp_vol7 = generate_vol_seq(70.0, 40.0);*/
	/*std::vector<double> temp_vol3 = generate_vol_seq(0.0, 130.0);
	std::vector<double> temp_vol4 = generate_vol_seq(0.0, 110.0);
	std::vector<double> temp_vol5 = generate_vol_seq(0.0, 80.0); */

	//temp_vol1 = temp_vol1 + temp_vol2 + temp_vol3 + temp_vol4 + temp_vol5 + temp_vol6 + temp_vol7;

	//const std::vector<double>  e1 = predict_disp(parameter_50, temp_vol1);


	std::vector<double> temp_disp1 = generate_disp_seq(0.0, 10.0);
	std::vector<double> temp_disp2 = generate_disp_seq(10.0, 0.0);
	std::vector<double> temp_disp3 = generate_disp_seq(0.0, 8.0);
	std::vector<double> temp_disp4 = generate_disp_seq(8.0, 5.0);

	temp_disp1 = temp_disp1 + temp_disp2 + temp_disp3 + temp_disp4;

	const std::vector<double> u = feed_forward(parameter_50, temp_disp1);
	const std::vector<double>  e1 = predict_disp(parameter_50, u);

	/* 绘图 */
	//plot(temp_vol1);
	plot(temp_disp1, e1);
	//plot(temp_vol1, e1);
	//plot(temp_vol2, e2);
	//plot(temp_vol3, e3);
	//plot(temp_vol4, e4);
	show();
}