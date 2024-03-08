/*<<Методы решений ОДУ>>
 @author: Шаманов Иван
 @author: Климов Олег
 2024*/
#pragma once
#include<map>
#include"LinOp.h"
#include"FileIO.h"
#include"NLEmethods.h"
//TODO:явный метод Эйлера;
template<typename DT, typename F>
bool EulerSolve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "eulersimple-1")
{
	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: Starting EulerSolve" << endl;
	cout << "log[INFO]: Opening a file to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> y_i = u_0;
		/*fpoints << t_i << endl;
		writeVectorToFile(fpoints, y_i);*/
		int ind = 0;
		vector<DT> y_ipp = u_0;
		while (t_i <= T)
		{
			y_ipp = y_i + tau * func(t_i, y_i);
			fpoints << t_i << endl;
			writeVectorToFile(fpoints, y_i);
			y_i = y_ipp;
			t_i += tau;
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}
template<typename DT, typename F>
vector<DT> k1(F func, DT tn, vector<DT> yn)
{
	return func(tn, yn);
}
template<typename DT, typename F>
vector<DT> k2(F func, DT tn, DT tau, vector<DT> yn)
{
	return func(tn + tau / 2, yn + 0.5 * tau * k1(func, tn, yn));
}
template<typename DT, typename F>
vector<DT> k2(F func, DT tn, DT tau, vector<DT> yn, vector<DT> k1)
{
	return func(tn + tau / 2, yn + 0.5 * tau * k1);
}
//TODO:двухстадийный метод 2 порядка точности;
template <typename DT, typename F>
bool RungeKutta2Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "rungekutta2-1")
{
	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: RungeKutta2Solve" << endl;
	cout << "log[INFO]: Opening a file to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> y_i = u_0;
		/*fpoints << t_i << endl;
		writeVectorToFile(fpoints, y_i);*/
		int ind = 0;
		vector<DT> y_ipp = u_0;
		while (t_i <= T)
		{
			y_ipp = y_i + tau * k2(func, t_i, tau, y_i);
			fpoints << t_i << endl;
			writeVectorToFile(fpoints, y_i);
			y_i = y_ipp;
			t_i += tau;
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}
template<typename DT, typename F>
vector<DT> k3(F func, DT tn, DT tau, vector<DT> yn)
{
	return func(tn + tau / 2, yn + 0.5 * tau * k2(func, tn, tau, yn));
}
template<typename DT, typename F>
vector<DT> k3(F func, DT tn, DT tau, vector<DT> yn, vector<DT> k2)
{
	return func(tn + tau / 2, yn + 0.5 * tau * k2);
}
template<typename DT, typename F>
vector<DT> k4(F func, DT tn, DT tau, vector<DT> yn)
{
	return func(tn + tau, yn + tau * k3(func, tn, tau, yn));
}
template<typename DT, typename F>
vector<DT> k4(F func, DT tn, DT tau, vector<DT> yn, vector<DT> k3)
{
	return func(tn + tau, yn + tau * k3);
}
//TODO:четырёхстадийный метод 4 порядка точности;
template <typename DT, typename F>
bool RungeKutta4Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "rungekutta4-1")
{
	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: RungeKutta4Solve" << endl;
	cout << "log[INFO]: Opening a file to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> y_i = u_0;/*
		fpoints << t_i << endl;
		writeVectorToFile(fpoints, y_i);*/
		int ind = 0;
		vector<DT> y_ipp = u_0;
		while (t_i <= T)
		{
			vector<DT> K1 = k1(func, t_i, y_i);
			vector<DT> K2 = k2(func, t_i, tau, y_i, K1);
			vector<DT> K3 = k3(func, t_i, tau, y_i, K2);
			vector<DT> K4 = k4(func, t_i, tau, y_i, K3);
			y_ipp = y_i + (tau/6) * (K1 + 2*K2 + 2*K3 + K4);
			fpoints << t_i << endl;
			writeVectorToFile(fpoints, y_i);
			y_i = y_ipp;
			t_i += tau;
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}
//TODO:неявный метод Эйлера;
template <typename DT, typename F>
bool ImplicitEulerSolve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "eulerImplicit-1")
{
	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: Starting ImplicitEulerSolve" << endl;
	cout << "log[INFO]: Opening a file to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> range{-10000000., 10000000.};
		vector<DT> y_i = u_0;
		/*fpoints << t_i << endl;
		writeVectorToFile(fpoints, y_i);*/
		int ind = 0;
		vector<DT> y_ipp = u_0;
		while (t_i <= T)
		{
			auto f1 = [&](vector<DT> y_ip) {return (1/ tau)*(y_ip-y_i) - func(t_i+tau, y_ip);};
			y_ipp = NewtonSolve(f1, range, 1e-6, 1000, y_i, false);
			fpoints << t_i << endl;
			writeVectorToFile(fpoints, y_i);
			y_i = y_ipp;
			t_i += tau;
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}
//TODO:симметричная схема;
template<typename DT, typename F>
bool SimSchemeSolve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "simScheme-1")
{
	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: Starting SimSchemeSolve" << endl;
	cout << "log[INFO]: Opening a file to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> range{ -10000000., 10000000. };
		vector<DT> y_i = u_0;
		/*fpoints << t_i << endl;
		writeVectorToFile(fpoints, y_i);*/
		int ind = 0;
		vector<DT> y_ipp = u_0;
		while (t_i <= T)
		{
			auto f1 = [&](vector<DT> y_ip) {return (1 / tau) * (y_ip - y_i) - 0.5 * (func(t_i, y_i) + func(t_i+tau, y_ip));};
			y_ipp = NewtonSolve(f1, range, 1e-6, 1000, y_i, false);
			fpoints << t_i << endl;
			writeVectorToFile(fpoints, y_i);
			y_i = y_ipp;
			t_i += tau;
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}
//TODO:явный четырёхшаговый метод Адамса—Башфорта 4 порядка точности;
template<typename DT,  typename F>
bool Adams4Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "adams4-1")
{
	string path = "OutputData\\" + filename;
	ofstream fpoints_init(path);
	cout << "log[INFO]: Adams4Solve" << endl;
	cout << "log[INFO]: Opening a file to write..." << endl;
	if (fpoints_init.is_open())
	{
		fpoints_init.close();
		RungeKutta4Solve(func, t_0, t_0+3*tau, tau, u_0, filename);
		fstream fpoints;
		fpoints.open(path, std::ios::app);
		string tmp_line;
		DT t_i = t_0;
		vector<vector<DT>> init_vals(readInitVals<DT>(path));
		vector<DT> y_m3(init_vals[1]);
		vector<DT> f_m3 = func(t_i, y_m3);
		vector<DT> y_m2(init_vals[3]);
		t_i += tau;
		vector<DT> f_m2 = func(t_i, y_m2);
		vector<DT> y_m1(init_vals[5]);
		t_i += tau;
		vector<DT> f_m1 = func(t_i, y_m1);
		vector<DT> y_i(init_vals[7]);
		t_i += tau;
		vector<DT> y_ipp = y_i;
		vector<DT> f_i = f_m1;
		while (t_i < T)
		{
			f_i = func(t_i, y_i);
			y_ipp = y_i + (tau / 24) * (55*f_i + (-59) * f_m1 + 37 * f_m2 +(-9) * f_m3);
			f_m3 = f_m2;
			f_m2 = f_m1;
			f_m1 = f_i;
			y_i = y_ipp;
			t_i += tau;
			fpoints << t_i << endl;
			writeVectorToFile(fpoints, y_ipp);
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}
//TODO:метод предиктор–корректор 4 порядка точности;
template<typename DT, typename F>
bool PredictCorrect4Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "predict4-1")
{
	string path = "OutputData\\" + filename;
	ofstream fpoints_init(path);
	cout << "log[INFO]: PredictCorrect4Solve" << endl;
	cout << "log[INFO]: Opening a file to write..." << endl;
	if (fpoints_init.is_open())
	{
		fpoints_init.close();
		RungeKutta4Solve(func, t_0, t_0 + 3 * tau, tau, u_0, filename);
		fstream fpoints;
		fpoints.open(path, std::ios::app);
		string tmp_line;
		DT t_i = t_0;
		vector<vector<DT>> init_vals(readInitVals<DT>(path));
		vector<DT> y_m3(init_vals[1]);
		vector<DT> f_m3 = func(t_i, y_m3);
		vector<DT> y_m2(init_vals[3]);
		t_i += tau;
		vector<DT> f_m2 = func(t_i, y_m2);
		vector<DT> y_m1(init_vals[5]);
		t_i += tau;
		vector<DT> f_m1 = func(t_i, y_m1);
		vector<DT> y_i(init_vals[7]);
		t_i += tau;
		vector<DT> y_ipp_0 = y_i;
		vector<DT> y_ipp = y_ipp_0;
		vector<DT> f_ipp = f_m1;
		vector<DT> f_i = f_m1;
		while (t_i < T)
		{
			f_i = func(t_i, y_i);
			y_ipp_0 = y_i + (tau / 24) * (55 * f_i + (-59) * f_m1 + 37 * f_m2 + (-9) * f_m3);
			t_i += tau;
			f_ipp = func(t_i, y_ipp_0);
			y_ipp = y_i + (tau / 24) * (9 * f_ipp + 19 * f_i + (-5)*f_m1 +  f_m2);
			f_m3 = f_m2;
			f_m2 = f_m1;
			f_m1 = f_i;
			y_i = y_ipp;
			fpoints << t_i << endl;
			writeVectorToFile(fpoints, y_ipp);
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}