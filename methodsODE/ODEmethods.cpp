#include "ODEmethods.h"


using namespace std;


/* Метод Эйлера Явный */
template<typename DT, typename F>
bool EulerSolve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename) {
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
			//fpoints << t_i << endl;
			writeVectorToFile(fpoints, t_i, y_i);
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


/* Метод Эйлера Неявный */
template <typename DT, typename F>
bool ImplicitEulerSolve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename) {
	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: Starting ImplicitEulerSolve" << endl;
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
			auto f1 = [&](vector<DT> y_ip) {return (1 / tau) * (y_ip - y_i) - func(t_i + tau, y_ip); };
			y_ipp = NewtonSolve(f1, range, 1e-6, 1000, y_i, false);
			//fpoints << t_i << endl;
			writeVectorToFile(fpoints, t_i, y_i);
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

/* Метод Симметричной схемы */
template<typename DT, typename F>
bool SimSchemeSolve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename) {
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
			auto f1 = [&](vector<DT> y_ip) {return (1 / tau) * (y_ip - y_i) - 0.5 * (func(t_i, y_i) + func(t_i + tau, y_ip)); };
			y_ipp = NewtonSolve(f1, range, 1e-6, 1000, y_i, false);
			//fpoints << t_i << endl;
			writeVectorToFile(fpoints, t_i, y_i);
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

/*
template<typename DT, typename F>
vector<DT> k1(F func, DT tn, vector<DT> yn) {
	return func(tn, yn);
}


template<typename DT, typename F>
vector<DT> k2(F func, DT tn, DT tau, vector<DT> yn) {
	return func(tn + tau / 2, yn + 0.5 * tau * k1(func, tn, yn));
}


template<typename DT, typename F>
vector<DT> k2(F func, DT tn, DT tau, vector<DT> yn, vector<DT> k1)
{
	return func(tn + tau / 2, yn + 0.5 * tau * k1);
}
*/

/* Метод Рунге-Кутты двустадийный */
template <typename DT, typename F>
bool RungeKutta2Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename) {

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
		vector<DT> k1 = u_0 , k2 = u_0;

		while (t_i <= T) {

			k1 = func(t_i, y_i);
			k2 = func(t_i + tau / 2, y_i + 0.5 * tau * k1);

			y_ipp = y_i + (0.5 * tau) * (k1 + k2);


			//fpoints << t_i << endl;
			writeVectorToFile(fpoints, t_i, y_i);
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


/* Метод Рунге-Кутты двустадийный с Автоматическим шагом */
template <typename DT, typename F>
bool RungeKutta2SolveAuto(F func, DT t_0, DT T, DT tau0, vector<DT> u_0, string filename) {


	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: RungeKutta2SolveAuto" << endl;
	cout << "log[INFO]: Opening a file to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> y_i = u_0;
		/*fpoints << t_i << endl;
		writeVectorToFile(fpoints, y_i);*/
		int ind = 0;
		vector<DT> y_ipp_1 = u_0, y_ipp_2 = u_0, y_ipp_3 = u_0;
		vector<DT> k1 = u_0, k2 = u_0;

		DT tau = tau0;

		while (t_i <= T) {

			/* Вычисляем следующую точку по половинному, целому и двойному шагу */

			// Вычисляем компоненты k при половинном шаге
			k1 = func(t_i + tau / 2, y_i);
			k2 = func(t_i + tau / 2, y_i + tau * k1);
			y_ipp_1 = y_i + (tau / 4.) * (k1 + k2);

			// Вычисляем компоненты k при целом шаге
			k1 = func(t_i + tau, y_i);
			k2 = func(t_i + tau, y_i + tau * k1);
			y_ipp_2 = y_i + (tau / 2.) * (k1 + k2);

			// Вычисляем компоненты k при двойном шаге
			/*
			k1 = func(t_i + 2 * tau, y_i);
			k2 = func(t_i + 2 * tau, y_i + 2 * tau * k1);
			y_ipp_3 = y_i + 2 * tau * (k1 + k2);
			*/
			
			//fpoints << t_i << endl;
			writeVectorToFile(fpoints, t_i, y_i);

			t_i += tau;

			if ((vec_norm((1. / 3) * (y_ipp_1 - y_ipp_2)) < tau0) /* and tau < tau0*/) {
				tau = 2 * tau;
			} else {
				tau = 0.5 * tau;
			}
			y_i = y_ipp_2;
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}


/* Метод Рунге-Кутты четырехстадийный*/
template <typename DT, typename F>
bool RungeKutta4Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename)
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
		vector<DT> k1 = u_0, k2 = u_0, k3 = u_0, k4 = u_0;
		while (t_i <= T) {

			k1 = func(t_i, y_i);
			k2 = func(t_i + tau / 2, y_i + 0.5 * tau * k1);
			k3 = func(t_i + tau / 2, y_i + 0.5 * tau * k2);
			k4 = func(t_i + tau, y_i + tau * k3);

			y_ipp = y_i + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
			//fpoints << t_i << endl;
			writeVectorToFile(fpoints, t_i, y_i);
			y_i = y_ipp;
			t_i += tau;
		}

		fpoints.close();
		return true;
	}

	else {
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	}
	return false;
}


/* Метод Рунге-Кутты четырехстадийный c автоматическим шагом */
template <typename DT, typename F>
bool RungeKutta4SolveAuto(F func, DT t_0, DT T, DT tau0, vector<DT> u_0, string filename) 
{
	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: RungeKutta4SolveAuto" << endl;
	cout << "log[INFO]: Opening a file to write..." << endl;
	if (fpoints.is_open())
	{
		DT time = t_0;
		vector<DT> y_i = u_0;/*
		fpoints << t_i << endl;
		writeVectorToFile(fpoints, y_i);*/
		int ind = 0;
		vector<DT> y_ipp_1 = u_0, y_ipp_2 = u_0, y_ipp_3 = u_0;
		vector<DT> k1 = u_0, k2 = u_0, k3 = u_0, k4 = u_0, K = u_0;

		DT eps_h = 1e-5;
		DT tau = tau0;
		while (time <= T) {

			// Вычисляем следующую точку по половинному, целому и двойному шагу

			// Вычисляем компоненты k при половинном шаге
			k1 = func(time + tau, y_i);
			k2 = func((time + tau / 4), y_i + (tau / 4) * k1);
			k3 = func((time + tau / 4), y_i + (tau / 4) * k2);
			k4 = func((time + tau / 2), y_i + (tau / 2) * k3);
			K = (1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
			y_ipp_1 = y_i + tau / 2 * K;

			// Вычисляем компоненты k при целом шаге
			k1 = func(time + tau, y_i);
			k2 = func((time + tau / 2), y_i + (tau / 2) * k1);
			k3 = func((time + tau / 2), y_i + (tau / 2) * k2);
			k4 = func((time + tau), y_i + tau * k3);
			K = (1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
			y_ipp_2 = y_i + tau * K;

			// Вычисляем компоненты k при двойном шаге
			/*
			k1 = func(time + 2 * tau, y_i);
			k2 = func((time + tau), y_i + tau * k1);
			k3 = func((time + tau), y_i + tau * k2);
			k4 = func((time + tau * 2), y_i + 2 * tau * k3);
			K = (1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
			y_ipp_3 = y_i + 2 * tau * K;
			*/
			time += tau;

			// Выбираем следующий используемый шаг
			if (vec_norm((1. / 15) * (y_ipp_1 - y_ipp_2)) < eps_h and tau < tau0 and tau > eps_h) {
				tau = 2 * tau;
			}

			else {
				tau = 0.5 * tau;
			}

			//fpoints << t_i << endl;
			y_i = y_ipp_2;
			writeVectorToFile(fpoints, time, y_i);
		
		
		}

		fpoints.close();
		return true;
	}

	else {
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	}
	return false;
}


/* Метод Адамса-Башфорта */
template<typename DT, typename F>
bool Adams4Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename) {

	string path = "OutputData\\" + filename;
	ofstream fpoints_init(path);
	cout << "log[INFO]: Adams4Solve" << endl;
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
		vector<DT> y_m3(init_vals[0]);
		y_m3.erase(y_m3.begin());
		vector<DT> f_m3 = func(t_i, y_m3);
		vector<DT> y_m2(init_vals[1]);
		y_m2.erase(y_m2.begin());
		t_i += tau;
		vector<DT> f_m2 = func(t_i, y_m2);
		vector<DT> y_m1(init_vals[2]);
		y_m1.erase(y_m1.begin());
		t_i += tau;
		vector<DT> f_m1 = func(t_i, y_m1);
		vector<DT> y_i(init_vals[3]);
		y_i.erase(y_i.begin());
		t_i += tau;
		vector<DT> y_ipp = y_i;
		vector<DT> f_i = f_m1;
		while (t_i < T)
		{
			f_i = func(t_i, y_i);
			y_ipp = y_i + (tau / 24) * (55 * f_i + (-59) * f_m1 + 37 * f_m2 + (-9) * f_m3);
			f_m3 = f_m2;
			f_m2 = f_m1;
			f_m1 = f_i;
			y_i = y_ipp;
			t_i += tau;
			//fpoints << t_i << endl;
			writeVectorToFile(fpoints, t_i, y_ipp);
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}


/* Метор Прогноз-Коррекция */
template<typename DT, typename F>
bool PredictCorrect4Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename) {

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
		vector<DT> y_m3(init_vals[0]);
		y_m3.erase(y_m3.begin());
		vector<DT> f_m3 = func(t_i, y_m3);
		vector<DT> y_m2(init_vals[1]);
		y_m2.erase(y_m2.begin());
		t_i += tau;
		vector<DT> f_m2 = func(t_i, y_m2);
		vector<DT> y_m1(init_vals[2]);
		y_m1.erase(y_m1.begin());
		t_i += tau;
		vector<DT> f_m1 = func(t_i, y_m1);
		vector<DT> y_i(init_vals[3]);
		y_i.erase(y_i.begin());
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
			y_ipp = y_i + (tau / 24) * (9 * f_ipp + 19 * f_i + (-5) * f_m1 + f_m2);
			f_m3 = f_m2;
			f_m2 = f_m1;
			f_m1 = f_i;
			y_i = y_ipp;
			//fpoints << t_i << endl;
			writeVectorToFile(fpoints, t_i, y_ipp);
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}
