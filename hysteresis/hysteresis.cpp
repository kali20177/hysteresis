#include "pch.h"
#include <Windows.h>
#include <matplotlibcpp.h>

using namespace matplotlibcpp;

#define M_PI 3.14159265358979323846

double interval = 0.0125;

double shape_control_fun(int data_num, int data_size, double* arr);

template <typename T>
std::vector<T>& operator+(std::vector<T>& v1, std::vector<T>& v2)
{
	v1.insert(v1.end(), v2.begin(), v2.end());
	return v1;
}

//产生两点之间序列
double temp = 150.0;

int count = 400;

std::vector<double> generate_vol_seq(double new_vol)
{
	std::vector<double> vol_seq;

	for (int i = 0; i < count; ++i)
	{
		double vol;
		if (new_vol >= temp)
			vol = fabs((new_vol - temp) / 2) * sin(2 * M_PI * i / (2 * count) - M_PI / 2) + fabs((new_vol - temp) / 2) + min(new_vol, temp);
		else
			vol = fabs((new_vol - temp) / 2) * sin(2 * M_PI * i / (2 * count) + M_PI / 2) + fabs((new_vol - temp) / 2) + min(new_vol, temp);

		vol_seq.emplace_back(vol);
	}
	temp = new_vol;
	return vol_seq;
}

double shape_control_fun(int data_num, int data_size, double* arr)
{
	double mid_pra = 0;
	double part = (double)data_num / data_size;

	if (part < 0.125)
	{
		mid_pra = *(arr + 3);
	}
	else if (part < 0.25)
	{
		mid_pra = *(arr + 4);
	}
	else if (part < 0.375)
	{
		mid_pra = *(arr + 5);
	}
	else if (part < 0.5)
	{
		mid_pra = *(arr + 6);
	}
	else if (part < 0.625)
	{
		mid_pra = *(arr + 7);
	}
	else if (part < 0.75)
	{
		mid_pra = *(arr + 8);
	}
	else if (part < 0.875)
	{
		mid_pra = *(arr + 9);
	}
	else
	{
		mid_pra = *(arr + 10);
	}

	return mid_pra;
}

std::vector<double> predict_disp(double* parameter, const std::vector<double>& input_vol)
{
	double k, a, y, eta = 0;
	double mid_prameter = 0;

	k = *parameter;
	a = *(parameter + 1);
	y = *(parameter + 2);
	eta = *(parameter + 11);

	std::vector<double> d_input;

	for (int i = 0; i < input_vol.size(); ++i)
	{
		if (i < input_vol.size() - 1)
			d_input.emplace_back((input_vol[i + 1] - input_vol[i]) / interval);
		else
			d_input.emplace_back(d_input[i - 1]);
	}

	return d_input;
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

	std::vector<double> temp_vol1 = generate_vol_seq(40.0);
	std::vector<double> temp_vol2 = generate_vol_seq(60.0);
	std::vector<double> temp_vol3 = generate_vol_seq(100.0);
	std::vector<double> temp_vol4 = generate_vol_seq(120.0);
	std::vector<double> temp_vol5 = generate_vol_seq(80.0);

	temp_vol1 = temp_vol1 + temp_vol2 + temp_vol3 + temp_vol4 + temp_vol5;

	const std::vector<double>  d_input = predict_disp(parameter_50, temp_vol1);
	plot(temp_vol1);
	plot(d_input);
	show();
}