#include <iostream>
#include<vector>
#include <iomanip>
#include<cmath>
#include <random>
#include<string>
#include<map>
double fun_error(double& x, double& c, double& d, double& error)
{
	return c * x + d + error;
}
double fun(double& x, double& c, double& d)
{
	return c * x + d;
}
double MNK(std::vector<std::pair<double, double>>& pair, double& N)
{
	double sum_1 = 0, sum_x = 0, sum_t = 0, sum_x2 = 0, sum2x = 0;

	for (size_t i = 0; i < N + 1; i++)
	{
		sum_1 += pair[i].first * pair[i].second;
		sum_x += pair[i].first;
		sum_t += pair[i].second;
		sum_x2 += pow(pair[i].first, 2);
	}

	return (N * sum_1 - sum_x * sum_t) / (N * sum_x2 - pow(sum_x, 2));
}
double passive_search(std::vector<std::pair<double, double>> & pair, double& N, double& c_d)
{
	double max_d = 2.5, min_d = -0.5, number = 1000;
	std::vector<double> vec_d;
	std::vector<double> sum;
	for (size_t i = 0; i < number; i++)
	{
		double sum_1 = 0;
		vec_d.push_back(min_d + ((max_d - min_d) / (number - 1)) * (i - 1));
		for (auto j = 0; j < N + 1; j++) {

			double new_y = fun(pair[j].first, c_d, vec_d[i]);
			sum_1 += pow(new_y - pair[j].second, 2);
		}
		sum.push_back(sum_1);

	}
	auto min_sum = std::min_element(sum.begin(), sum.end());
	auto num = std::distance(sum.begin(), min_sum);
	return vec_d[num];
}

int main()
{
	double c = 3, d = 1, a = -1, b = 3, N = 16, A = 3, x;
	std::vector<std::pair<double, double>> pair;
	double step = double((b - a)) / N;
	std::cout << step << std::endl;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> error(-0.5, 0.5);
	double my_rd;
	for (size_t i = 0; i < N + 1; i++)
	{
		x = a + i * step;
		double e = error(gen);
		std::pair<double, double> x_y_step;
		x_y_step.first = x;
		x_y_step.second = fun_error(x, c, d, e);
		pair.push_back(x_y_step);
		std::cout << x_y_step.first << "  " << x_y_step.second << std::endl;

	}
	double c_d = 0;
	c_d = MNK(pair, N);
	std::cout << "-----------------------------------" << std::endl;
	std::cout << "c=  " << c_d << std::endl;
	double d_d = 0;
	d_d = passive_search(pair, N, c_d);
	std::cout << "d=  " << d_d << std::endl;
	for (size_t i = 0; i < N + 1; i++)
	{
		x = a + i * step;
		std::pair<double, double> x_y_step;
		x_y_step.first = x;
		x_y_step.second = fun(x, c_d, d_d);
		pair.push_back(x_y_step);
		std::cout << x_y_step.first << "  " << x_y_step.second << std::endl;
	}
}