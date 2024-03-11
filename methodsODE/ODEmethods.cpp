#include "ODEmethods.h"

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
