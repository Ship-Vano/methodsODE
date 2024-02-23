/*<<Лабораторная работа №1. Методы решения ОДУ>>
 @author: Шаманов Иван
 2024*/
#include<iostream>
#include"ODEmethods.h"
int main()
{
	auto fn1 = [](double x, vector<double> u) {return vector<double>{2*u[0] + u[1]*u[1] -1, 
																	6*u[0]-u[1]*u[1]+1};};
	double t_0 = 0;
	double T = 2;
	double tau = 0.001;
	vector<double> u_0{0.,0.};
	EulerSolve(fn1, t_0, T, tau, u_0);
	ImplicitEulerSolve(fn1, t_0, T, tau, u_0);
	RungeKutta2Solve(fn1, t_0, T, tau, u_0);
	RungeKutta4Solve(fn1, t_0, T, tau, u_0);
	return 0;
}