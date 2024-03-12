﻿#include "ODEmethods.h"

/* Метод Эйлера Явный */
template<typename DT, typename F>
bool EulerSolve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename) {
	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: Starting EulerSolve" << endl;
	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> y_i = u_0;
		/*fpoints << t_i << endl;
		writeVectorToFile(fpoints, y_i);*/
		int ind = 0;
		vector<DT> y_ipp = u_0;
		writeVectorToFile(fpoints, t_i, y_i);
		while (abs(T-t_i) >=  1e-8)
		{
			y_ipp = y_i + tau * func(t_i, y_i);
			//fpoints << t_i << endl;
			y_i = y_ipp;
			t_i += tau;
			writeVectorToFile(fpoints, t_i, y_i);
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
	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> range{ -10000000., 10000000. };
		vector<DT> y_i = u_0;
		int ind = 0;
		vector<DT> y_ipp = u_0;
		writeVectorToFile(fpoints, t_i, y_i);
		while (abs(T - t_i) >= 1e-8)
		{
			auto f1 = [&](vector<DT> y_ip) {return (1 / tau) * (y_ip - y_i) - func(t_i + tau, y_ip); };
			y_ipp = NewtonSolve(f1, range, 1e-6, 1000, y_i, false);

			y_i = y_ipp;
			t_i += tau;
			writeVectorToFile(fpoints, t_i, y_i);
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
	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> range{ -10000000., 10000000. };
		vector<DT> y_i = u_0;
		int ind = 0;
		vector<DT> y_ipp = u_0;
		writeVectorToFile(fpoints, t_i, y_i);
		while (abs(T - t_i) >= 1e-8)
		{
			auto f1 = [&](vector<DT> y_ip) {return (1 / tau) * (y_ip - y_i) - 0.5 * (func(t_i, y_i) + func(t_i + tau, y_ip)); };
			y_ipp = NewtonSolve(f1, range, 1e-6, 1000, y_i, false);
			
			y_i = y_ipp;
			t_i += tau;
			writeVectorToFile(fpoints, t_i, y_i);
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
	cout << "log[INFO]: Starting RungeKutta2Solve" << endl;
	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> y_i = u_0;
		vector<DT> y_ipp = u_0;
		vector<DT> k1 = u_0 , k2 = u_0;
		writeVectorToFile(fpoints, t_i, y_i);
		while (abs(T - t_i) > 1e-8) {
			//из методички:
			k1 = func(t_i, y_i);
			k2 = func(t_i + tau, y_i +  tau * k1);

			y_ipp = y_i + (0.5 * tau) * (k1 + k2);
			//из учебника:
			/*
			* k1 = func(t_i, y_i);
			* k2 = func(t_i + tau/2, y_i +  tau/2 * k1);
			* y_ipp = y_i + tau * k2;
			*/
			y_i = y_ipp;
			t_i += tau;
			writeVectorToFile(fpoints, t_i, y_i);			
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}


/* Метод Рунге-Кутты двустадийный с Автоматическим шагом */
template<typename Ff, typename DT>
DT dif_eval2(Ff func, DT &tau, DT t_i,  DT t_edge,  vector<DT> &k1, vector<DT> &k2, vector<DT> y_i, vector<DT> & y_ipp_1, vector<DT>& y_ipp_2)
{
	// Вычисляем компоненты k при половинном шаге 
	//первый половинный шаг
	DT tau_half = tau / 2;
	DT temp_t = t_i;
	vector<DT> y_ipp_1_0 = y_i;
	while (t_edge - temp_t > 0)
	{
		temp_t += tau_half;
		k1 = func(temp_t - tau_half, y_ipp_1_0);
		k2 = func(temp_t + tau_half, y_ipp_1_0 + tau_half * k1);
		y_ipp_1_0 = y_ipp_1_0 + (tau_half / 2) * (k1 + k2);
	}
	//cout << t_edge << "v/s" << temp_t << endl;
	//cout << "----t_edge = " << temp_t;
	y_ipp_1 = y_ipp_1_0;
	// Вычисляем компоненты k при целом шаге
	temp_t = t_i;
	y_ipp_1_0 = y_i;
	while (t_edge - temp_t > 0)
	{
		temp_t += tau;
		k1 = func(temp_t-tau, y_ipp_1_0);
		k2 = func(temp_t + tau, y_ipp_1_0 + tau * k1);
		y_ipp_1_0 = y_ipp_1_0 + (tau / 2) * (k1 + k2);
	}
	y_ipp_2 = y_ipp_1_0;
	//проверяем апостериорную погрешность
	int p = 2;
	DT denom = pow(2, p) - 1;
	DT difference = vec_norm(y_ipp_2 - y_ipp_1) / denom;
	//cout << "counted diff" << endl;
	return difference;
}
template <typename DT, typename F>
bool RungeKutta2SolveAuto(F func, DT t_0, DT T, DT tau0, vector<DT> u_0, string filename, DT eps = 1e-6) {


	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: Starting RungeKutta2SolveAuto" << endl;
	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> y_i = u_0;
		vector<DT> y_ipp_1_0 = u_0, y_ipp_1 = u_0, y_ipp_2 = u_0;
		vector<DT> k1 = u_0, k2 = u_0;

		DT tau = tau0;
		writeVectorToFile(fpoints, t_i, y_i);
		while (T - t_i > 0) {

			DT edge = t_i + tau;
			DT difference = dif_eval2(func, tau, t_i, edge, k1, k2, y_i, y_ipp_1, y_ipp_2);
			while (difference >= eps)
			{
				tau /= 2;
				difference = dif_eval2(func, tau, t_i, edge, k1, k2, y_i, y_ipp_1, y_ipp_2);
			}
			
			y_i = y_ipp_2;
			t_i += tau;
			writeVectorToFile(fpoints, t_i, y_i);
			if (difference <= eps * 1e-3)
			{
				tau *= 2;
			}
			continue;
			
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}

//template <typename DT, typename F>
//bool RungeKutta2SolveAuto(F func, DT t_0, DT T, DT tau0, vector<DT> u_0, string filename, DT eps = 1e-6) {
//
//	string path = "OutputData\\" + filename;
//	ofstream fpoints(path);
//	cout << "log[INFO]: Starting RungeKutta2SolveAuto" << endl;
//	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
//
//	if (fpoints.is_open()) {
//
//		DT t_i = t_0;
//		vector<DT> y_i = u_0;
//		vector<DT> y_ipp = u_0, y_ipp_1 = u_0, y_ipp_2 = u_0;
//		vector<DT> k1 = u_0, k2 = u_0;
//
//		DT tau = tau0; // Динамический шаг
//		DT err = eps; // Погрешность 
//		
//		// Запись первой точки (н.у.)
//		writeVectorToFile(fpoints, t_i, y_i);
//
//		// Стадийный процесс
//		while (T - t_i > 0) {
//			
//			// Вычисляем решение при обычном шаге
//			k1 = func(t_i, y_i);
//			k2 = func(t_i + tau, y_i + tau * k1);
//			y_ipp_1 = y_i + (0.5 * tau) * (k1 + k2);
//
//			// Вычисляем при половинном шаге (2 стадии)
//			k1 = func(t_i, y_i);
//			k2 = func(t_i + 0.5 * tau, y_i + 0.5 * tau * k1);
//			y_ipp_2 = y_i + (0.25 * tau) * (k1 + k2);
//			
//			k1 = func(t_i, y_ipp_2);
//			k2 = func(t_i + 0.5 * tau, y_ipp_2 + 0.5 * tau * k1);
//			y_ipp_2 = y_ipp_2 + (0.25 * tau) * (k1 + k2);
//
//			// Находим погрешность
//			err =  vec_norm_inf(y_ipp_2 - y_ipp_1) / 3;
//
//			// В зависимости от погрешности выбираем следующий шаг:
//			
//			
//		
//			//  1) Если слишком большая точность, то увеличиваем
//			if (err < eps * 1e-5) {
//
//				// Обновляем текущее решение
//				y_ipp = y_ipp_2;
//				y_i = y_ipp;
//				t_i += tau;
//
//				// Записываем в фалй
//				writeVectorToFile(fpoints, t_i, y_i);
//				continue;
//			}
//
//
//			// 2) Если погрешность больше точности, то уменьшаем шаг
//			
//			// Цикл: будем делать меньше шаг, пока не дойдем до нужной погрешности
//
//			if (err < eps) { // Бесполезный оператор, но пусть будет :)
//				while (err > eps) {
//
//					y_ipp_1 = y_ipp_2;
//					tau = tau / 2;
//
//					// Вычисляем при половинном шаге (2 стадии)
//					k1 = func(t_i, y_i);
//					k2 = func(t_i + 0.5 * tau, y_i + 0.5 * tau * k1);
//					y_ipp_2 = y_i + (0.25 * tau) * (k1 + k2);
//
//					k1 = func(t_i, y_ipp_2);
//					k2 = func(t_i + 0.5 * tau, y_ipp_2 + 0.5 * tau * k1);
//					y_ipp_2 = y_ipp_2 + (0.25 * tau) * (k1 + k2);
//
//					// Находим погрешность
//					err = vec_norm_inf(y_ipp_2 - y_ipp_1) / 3;
//				}
//
//				// Обновляем текущее решение
//				y_ipp = y_ipp_2;
//				y_i = y_ipp;
//				t_i += tau;
//
//				// Записываем в фалй
//				writeVectorToFile(fpoints, t_i, y_i);
//
//				continue;
//			}
//			else { 
//				// 3) Оставляем шаг каким он есть, если точность
//				// в диапазазоне между заданной и сверх
//
//				// Обновляем текущее решение
//				y_ipp = y_ipp_2;
//				y_i = y_ipp;
//				t_i += tau;
//
//				// Записываем в фалй
//				writeVectorToFile(fpoints, t_i, y_i);
//
//			}
//		}
//
//	} else
//		cout << "log[ERROR]: Couldn't open or create a file" << endl;
//	return false;
//
//}






/* Метод Рунге-Кутты четырехстадийный*/
template <typename DT, typename F>
bool RungeKutta4Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename)
{
	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: Starting RungeKutta4Solve" << endl;
	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> y_i = u_0;
		vector<DT> y_ipp = u_0;
		vector<DT> k1 = u_0, k2 = u_0, k3 = u_0, k4 = u_0;
		writeVectorToFile(fpoints, t_i, y_i);
		while (abs(T - t_i) > 1e-8) {

			k1 = func(t_i, y_i);
			k2 = func(t_i + tau / 2, y_i + 0.5 * tau * k1);
			k3 = func(t_i + tau / 2, y_i + 0.5 * tau * k2);
			k4 = func(t_i + tau, y_i + tau * k3);

			y_ipp = y_i + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

			y_i = y_ipp;
			t_i += tau;
			writeVectorToFile(fpoints, t_i, y_i);
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
template<typename Ff, typename DT>
DT dif_eval4(Ff func, DT& tau, DT t_i, DT t_edge, vector<DT>& k1, vector<DT>& k2, vector<DT>& k3, vector<DT>& k4, vector<DT> y_i, vector<DT>& y_ipp_1, vector<DT>& y_ipp_2)
{
	// Вычисляем компоненты k при половинном шаге 
	//первый половинный шаг
	DT tau_half = tau / 2;
	DT temp_t = t_i;
	vector<DT> y_ipp_1_0 = y_i;
	while (t_edge - temp_t > 0)
	{
		temp_t += tau_half;
		k1 = func(temp_t - tau_half, y_ipp_1_0);
		k2 = func(temp_t - tau_half / 2, y_ipp_1_0 + 0.5 * tau_half * k1);
		k3 = func(temp_t - tau_half / 2, y_ipp_1_0 + 0.5 * tau_half * k2);
		k4 = func(temp_t, y_ipp_1_0 + tau_half * k3);
		y_ipp_1_0 = y_ipp_1_0 + (tau_half / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

	}
	//cout << t_edge << "v/s" << temp_t << endl;
	//cout << "----t_edge = " << temp_t;
	y_ipp_1 = y_ipp_1_0;
	// Вычисляем компоненты k при целом шаге
	temp_t = t_i;
	y_ipp_1_0 = y_i;
	while (t_edge - temp_t > 0)
	{
		temp_t += tau;
		k1 = func(temp_t - tau, y_ipp_1_0);
		k2 = func(temp_t - tau / 2, y_ipp_1_0 + 0.5 * tau * k1);
		k3 = func(temp_t - tau / 2, y_ipp_1_0 + 0.5 * tau * k2);
		k4 = func(temp_t, y_ipp_1_0 + tau * k3);
		y_ipp_1_0 = y_ipp_1_0 + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
	}
	y_ipp_2 = y_ipp_1_0;
	//проверяем апостериорную погрешность
	int p = 4;
	DT denom = pow(2, p) - 1;
	DT difference = vec_norm(y_ipp_2 - y_ipp_1) / denom;
	//cout << "counted diff" << endl;
	return difference;
}
template <typename DT, typename F>
bool RungeKutta4SolveAuto(F func, DT t_0, DT T, DT tau0, vector<DT> u_0, string filename, DT eps=1e-6) 
{
	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	cout << "log[INFO]: Starting RungeKutta4SolveAuto" << endl;
	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
	if (fpoints.is_open())
	{
		DT t_i = t_0;
		vector<DT> y_i = u_0;
		int ind = 0;
		vector<DT> y_ipp_1 = u_0, y_ipp_2 = u_0;
		vector<DT> k1 = u_0, k2 = u_0, k3 = u_0, k4 = u_0, K = u_0;
		DT tau = tau0;
		writeVectorToFile(fpoints, t_i, y_i);
		while (T - t_i > 0) {

			DT edge = t_i + tau;
			DT difference = dif_eval4(func, tau, t_i, edge, k1, k2,k3, k4, y_i, y_ipp_1, y_ipp_2);
			while (difference >= eps)
			{
				tau /= 2;
				difference = dif_eval4(func, tau, t_i, edge, k1, k2, k3, k4, y_i, y_ipp_1, y_ipp_2);
			}
			if (difference < eps)
			{
				y_i = y_ipp_2;
				t_i += tau;
				writeVectorToFile(fpoints, t_i, y_i);
				if (difference <= eps * 1e-3)
				{
					tau *= 2;
				}
				continue;
			}
		
		
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
	cout << "log[INFO]: Starting Adams4Solve" << endl;
	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
	if (fpoints_init.is_open())
	{
		fpoints_init.close();
		DT rk_edge = t_0 + 3 * tau;
		if (T - rk_edge < 0)
		{
			cout << "log[ERROR]: Invalid net. Solving with \"RungeKutta4Solve\"" << endl;
			RungeKutta4Solve(func, t_0, T, tau, u_0, filename);
			return false;
		}
		RungeKutta4Solve(func, t_0, rk_edge, tau, u_0, filename);
		fstream fpoints;
		fpoints.open(path, std::ios::app);
		string tmp_line;
		DT t_i = t_0;
		vector<vector<DT>> init_vals(readInitVals<DT>(path));
		int init_size = init_vals.size();
		if (init_size < 4 || init_size > 4)
		{
			cout << "log[ERROR]: Invalid amount of start values" << endl;
			return false;
		}
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
		while (abs(T - t_i) > 1e-8)
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
	cout << "log[INFO]: Starting PredictCorrect4Solve" << endl;
	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
	if (fpoints_init.is_open())
	{
		fpoints_init.close();
		DT rk_edge = t_0 + 3 * tau;
		if (T - rk_edge < 0)
		{
			cout << "log[ERROR]: Invalid net. Solving with \"RungeKutta4Solve\"" << endl;
			RungeKutta4Solve(func, t_0, T, tau, u_0, filename);
			return false;
		}
		RungeKutta4Solve(func, t_0, rk_edge, tau, u_0, filename);
		fstream fpoints;
		fpoints.open(path, std::ios::app);
		string tmp_line;
		DT t_i = t_0;
		vector<vector<DT>> init_vals(readInitVals<DT>(path));
		int init_size = init_vals.size();
		if (init_size < 4 || init_size > 4)
		{
			cout << "log[ERROR]: Invalid amount of start values" << endl;
			return false;
		}
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
		while (abs(T - t_i) > 1e-8)
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

			writeVectorToFile(fpoints, t_i, y_ipp);
		}
		fpoints.close();
		return true;
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	return false;
}



/* Функции для постройки фазовой плоскостей */

// Метод, возвращающий массив
template <typename DT, typename F>
vector<vector<DT>> MethodRungeKutta4(F func, DT t_0, DT T, DT tau, vector<DT> u_0) {
	
	DT t_i = t_0;
	vector<DT> y_i = u_0;
	vector<DT> y_ipp = u_0;
	vector<DT> k1 = u_0, k2 = u_0, k3 = u_0, k4 = u_0;
	vector<vector<DT>> solve;
	solve.push_back(u_0);
	//cout << solve[0][0] << " " << solve[0][1] << endl;

	//int n = (T - t_0) / tau;
	
	for (double t_i = t_0 + tau; t_i <= T; t_i += tau) {

		/* Метод Рунге-Кутта 4-стадийный */
		//k1 = func(t_i, y_i);
		//k2 = func(t_i + tau / 2, y_i + 0.5 * tau * k1);
		//k3 = func(t_i + tau / 2, y_i + 0.5 * tau * k2);
		//k4 = func(t_i + tau, y_i + tau * k3);
		//y_ipp = y_i + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

		/* Метод Эйлера явный */
		y_ipp = y_i + tau * func(t_i, y_i);

		y_i = y_ipp;
	
		solve.push_back(y_ipp);

	}

	return solve;
}


template <typename DT>
bool build_faze() {

	ofstream file("OutputData\\FZ1");
	if (file.is_open() == false) {
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
		return false;
	}

	cout << "log[INFO]: Start make fazoviy face :) " << endl;
	// Уравнение, которое решаем
	auto fn = [](double x, vector<double> u) {
		return vector<double>{2 * u[0] + u[1] * u[1] - 1, 6 * u[0] - u[1] * u[1] + 1};
	};

	// Задача из 16 варианта
	auto fn16 = [](double x, vector<double> u) {

		double B = 0.16;
		double alpha = 0.05;
		double gamma = 0.5;
		double w = 0.833;
		double k = -0.5;

		return vector<double>{u[1], B * cos(w * x) - alpha * u[1] - k * u[0] - gamma * u[0] * u[0] * u[0]};
		};

	// Область
	vector<vector<DT>> diapazon = { {-2, 2}, { -2, 2} };


	/* График сильно менятеся в зависимости от следующих параметров */

	// Шаг по области
	DT h = 0.1;

	// Шаг в решении
	DT tau = 0.1;

	// Диапазон времени в решении
	vector<DT> time_diapazon = { 0., 0.9 };



	// Решение в точке (x, y)
	vector<vector<DT>> local_solve;

	for (double x = diapazon[0][0]; x <= diapazon[0][1]; x += h) {
		for (double y = diapazon[1][0]; y <= diapazon[1][1]; y += h) {

			// Решаем в точке (x, y)
			local_solve = MethodRungeKutta4(fn16, time_diapazon[0], time_diapazon[1], tau, { x, y });
			//file << local_solve[0][0];

			// Запись U1
			for (int i = 0; i < local_solve.size(); i++) {
				file << local_solve[i][0] << " ";
			}
			file << endl;

			 //Запись U2
			for (int i = 0; i < local_solve.size(); i++) {
				file << local_solve[i][1] << " ";
			}
			file << endl;

		}
	}
	file.close();
	return true;
}