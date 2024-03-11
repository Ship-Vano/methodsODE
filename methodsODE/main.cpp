/*<<Лабораторная работа №1. Методы решения ОДУ>>
 @author: Шаманов Иван
 @author: Климов Олег
 2024*/

#include <iostream>
#include "ODEmethods.cpp"
using namespace std;


/* Функция выполнения всех методов */
template<typename F>
void experiment(F func, string filename, string fprefix, double q = 0.5)
{
	cout << "-----------" << endl;
	cout << "log[INFO]: Conducting an experiment on " << fprefix << endl;
	cout << "-----------" << endl;
	vector<vector<double>> init_vals(readInitVals<double>(filename));
	int init_size = init_vals.size();
	if (init_size < 4)
	{
		cout << "log[ERROR] Invalid InputData file (init file must include all start values)" << endl;
	}
	double t_0 = init_vals[0][0];
	cout << "log[INFO] Starting time (t_0) is set to " << t_0 << endl;
	double T = init_vals[1][0];
	cout << "log[INFO] Edge time (T) is set to " << T << endl;
	double tau = init_vals[2][0];
	cout << "log[INFO] Step (tau) is set to " << tau << endl;
	vector<double> u_0;
	for (int j = 0; j < init_vals[3].size(); ++j) {
		u_0.push_back(init_vals[3][j]);
	}
	cout << "log[INFO] Start point (u_0) is set to: "; out(u_0);
	double tolerance = 1e-2;
	if (init_size > 4)
	{
		tolerance = init_vals[4][0];
	}
	cout << "log[INFO] Tolerance for auto-step methods is set to " << tolerance << endl;
	cout << "-----------" << endl;
	for (int i = 0; i < 5; ++i) {
		cout << "log[INFO] Evaluating on the net with tau = " << tau << endl;
		cout << "-----------" << endl;
		EulerSolve(func, t_0, T, tau, u_0, fprefix + "-eulersimple-" + to_string(i));
		ImplicitEulerSolve(func, t_0, T, tau, u_0, fprefix + "-eulerImplicit-" + to_string(i));
		RungeKutta2Solve(func, t_0, T, tau, u_0, fprefix + "-rungekutta2-" + to_string(i));
		RungeKutta2SolveAuto(func, t_0, T, tau, u_0, fprefix + "-rungekutta2auto-" + to_string(i), tolerance);
		RungeKutta4Solve(func, t_0, T, tau, u_0, fprefix + "-rungekutta4-" + to_string(i));
		RungeKutta4SolveAuto(func, t_0, T, tau, u_0, fprefix + "-rungekutta4auto-" + to_string(i), tolerance);
		SimSchemeSolve(func, t_0, T, tau, u_0, fprefix + "-simScheme-" + to_string(i));
		Adams4Solve(func, t_0, T, tau, u_0, fprefix + "-adams4-" + to_string(i));
		PredictCorrect4Solve(func, t_0, T, tau, u_0, fprefix + "-predict4-" + to_string(i));
		tau *= q;
		cout << "-----------" << endl;
	}

	cout << "log[INFO]: Completed the experiment on " << fprefix << " !" << endl;
}


int main() {

	// Тест 0 (Маятник)
	auto fn0 = [](double t, vector<double> u) { return vector<double>{ u[1], -1 * u[0]};};
	experiment(fn0, "InputData\\task0_1", "f0", 0.5);


	// Тест 1
	auto fn1 = [](double x, vector<double> u) { return vector<double>{2 * u[0] + u[1] * u[1] - 1, 6 * u[0] - u[1] * u[1] + 1}; };
	//experiment(fn1, "InputData\\task1_1", "f1", 0.5);
	
	// Тест 2
	auto fn2 = [](double x, vector<double> u) { return vector<double>{1 - u[0]*u[0] - u[1] * u[1], 2*u[0]};};
	//experiment(fn2, "InputData\\task2_1", "f2", 0.5);
	
	// Тест 3
	auto fn3 = [](double x, vector<double> u) { return vector<double>{10*(u[1]-u[0]), u[0]*(28-u[2])-u[1], u[0]*u[1] - 8 / 3 * u[2]};};
	//experiment(fn3, "InputData\\task3_1", "f3", 0.5);

	return 0;
}

//TODO: ввод начальных данных из файла +
//TODO: реализация тестовых функций (лямбды + эталонные решения) +-
//TODO: оценка погрешностей (5 сеток, таблица) +
//TODO: правило Эйткена -
//TODO: контроль шага в Р-К +
//TODO: 4 графика для стратегии контроля шага -
//TODO: маятник +
//TODO:	особые точки -
