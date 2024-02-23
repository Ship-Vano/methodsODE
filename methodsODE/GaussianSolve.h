#pragma once
#include"LinOp.h"
template<typename LT>
class System
{
public:
	System();
	System(vector<vector<LT>> A, vector<LT> b);
	System(string file_loc);
	void output();
	void sysoutput();
	void soloutput();
	/*vector<LT> error_vector();
	vector<LT> SimpleIterationSolve(vector<LT>& x_0, LT EPS);
	vector<LT> JacobySolve(vector<LT>& x_0, LT EPS);
	vector<LT> SeidelSolve(vector<LT>& x_0, LT EPS);
	vector<LT> RelaxSolve(vector<LT>& x_0, LT EPS);
	vector<LT> TripleBigRelaxSolve(vector<vector<LT>>& Sys, vector<LT>& x_0, LT EPS);*/
	int n = 0;
	vector<vector<LT>> MatrixA;
	vector<LT> RightVect;
	vector<LT> SolutionX;
	vector<LT> AnaliticalSol;
	vector<LT> X0;
	struct LogInfo
	{
		int iterations = 0;
		LT tau;
		LT C_norm;
		LT achieved_prec;
		LT error = 0;
		int aprior;
		int apost;
		vector<LT> x_0;
		vector<LT> error_vector;
	};

	LogInfo SimpIt_log_info;
	LogInfo Jacob_log_info;
	LogInfo Seidel_log_info;
	LogInfo Relax_log_info;
	bool GaussianPartChoiceSolve(LT EPS);
public:
	friend vector<vector <LT>> InverseMatrix(vector<vector <LT>> matrix, LT EPS);
protected:
	vector<vector<LT>> SystemMatrix;
	void file_init(string file_loc); // функция инициализации через файл
	int max_i(int index, LT EPS);
	void coef_reduce(int index, LT coef);
	void str_reduce(int current_ind, int upper_ind);


};

/*КОНСТРУКТОРЫ*/

template<typename LT>
System<LT>::System()  // конструктор из файла
{

}

template<typename LT>
System<LT>::System(vector<vector<LT>> A, vector<LT> b)  // конструктор из файла
{
	MatrixA = A;
	RightVect = b;
	n = MatrixA.size();
	for (int i = 0; i < n; i++)
	{
		SystemMatrix.push_back(MatrixA[i]);
		SystemMatrix[i].push_back(RightVect[i]);
	}

}

template<typename LT>
System<LT>::System(string file_loc)  // конструктор из файла
{
	file_init(file_loc);
}

template<typename LT>
void System<LT>::file_init(string file_loc)
{
	ifstream file(file_loc); // открываем файл
	if (file.is_open()) // вызов метода is_open()
	{
		string tmp_line;
		//file >> tmp_line;
		file >> tmp_line;

		cout << "Тест №" << tmp_line << endl;
		while (getline(file, tmp_line))
		{
			istringstream ss(tmp_line);
			SystemMatrix.emplace_back(istream_iterator<double>(ss), \
				istream_iterator<double>());
		}
		SystemMatrix.erase(SystemMatrix.begin());
		cout << "Успешно создали систему!!!" << endl;
		file.close();
		n = SystemMatrix[0].size() - 1;
		cout << "Размерность матрицы А: " << n << endl;
		//заполнение матрицы системы и вектора правой части:
		for (int i = 0; i < n; i++)
		{
			MatrixA.push_back({});
			for (int j = 0; j <= n; j++)
			{
				if (j == n)
				{
					RightVect.push_back(SystemMatrix[i][j]);
					continue;
				}
				MatrixA[i].push_back(SystemMatrix[i][j]);
			}
		}
	}
	else
	{
		cout << "Файл не открыт!\n\n" << endl;
	}
}

/*ВЫВОДЫ*/

template<typename LT>
void System<LT>::output()
{
	string sep1 = "_____________________________________________";
	string sep2 = "---------------------";
	cout << sep1 << endl;
	cout << "\n Система выглядит следующим образом: " << endl;
	cout << sep2 << endl;
	cout << "Матрица А: \n" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << fixed << setprecision(6) << setw(12) << setfill(' ') << MatrixA[i][j] << "  ";
		}
		cout << endl;
	}
	cout << sep2 << endl;

	cout << "Вектор правой части: \n" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << fixed << setprecision(6) << setw(12) << setfill(' ') << RightVect[i] << endl;
	}

	cout << sep1 << endl;
}

template<typename LT>
void System<LT>::sysoutput()
{
	cout << "Матрица системы: \n" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			cout << fixed << setprecision(6) << setw(12) << setfill(' ') << SystemMatrix[i][j] << "  ";
		}
		cout << endl;
	}
	cout << "-----------" << endl;
}

template<typename LT>
void System<LT>::soloutput()
{
	cout << "Вектор-решение Х: \n" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << scientific << setprecision(12) << setw(24) << setfill(' ') << SolutionX[i] << endl;
	}
}


template<typename LT>
int System<LT>::max_i(int index, LT EPS)
{
	int max_i = index;
	LT max_ai = fabs(SystemMatrix[index][index]);
	for (int i = index; i < n; i++)
		if (fabs(SystemMatrix[i][index]) > max_ai)
		{
			max_i = i;
			max_ai = fabs(SystemMatrix[i][index]);
		}
	//if (max_ai < EPS) return -1;
	return max_i;
}

template<typename LT>
void System<LT>::coef_reduce(int index, LT coef)
{
	for (int j = 0; j < n + 1; j++)
		SystemMatrix[index][j] /= coef;
}

template<typename LT>
void System<LT>::str_reduce(int current_ind, int upper_ind)
{
	LT coef = SystemMatrix[current_ind][upper_ind];
	for (int i = 0; i < n + 1; i++)
		SystemMatrix[current_ind][i] -= coef * SystemMatrix[upper_ind][i];
}

template<typename LT>
bool System<LT>::GaussianPartChoiceSolve(LT EPS)
{
	//ПРЯМОЙ ХОД:
	for (int i = 0; i < n; i++)
	{
		int max_pos = max_i(i, EPS);
		if (max_pos != -1)
		{
			SystemMatrix[i].swap(SystemMatrix[max_pos]);
			coef_reduce(i, SystemMatrix[i][i]);
			if (i == n) break;

			for (int j = i + 1; j < n; j++)
			{
				str_reduce(j, i);
			}

		}
		else {
			cout << "Матрица вырождена" << endl;
			return false;
		}
	}
	//ОБРАТНЫЙ ХОД:
	SolutionX.push_back(SystemMatrix[n - 1][n]);
	for (int i = n - 2; i > -1; i--)
	{
		LT xi = SystemMatrix[i][n];
		for (int j = i + 1; j < n; j++)
		{
			xi -= SystemMatrix[i][j] * SolutionX[n - j - 1];
		}
		SolutionX.push_back(xi);
	}
	reverse(SolutionX.begin(), SolutionX.end());
	return true;
}
template<typename LT>
vector<vector <LT>> InverseMatrix(vector<vector <LT>> matrix, LT EPS)
{
	int n = matrix.size();

	vector<LT> e1({ 1,0,0,0 });
	vector<LT> e2({ 0,1,0,0 });
	vector<LT> e3({ 0,0,1,0 });
	vector<LT> e4({ 0,0,0,1 });

	System<LT> A1(matrix, e1);
	System<LT> A2(matrix, e2);
	System<LT> A3(matrix, e3);
	System<LT> A4(matrix, e4);

	A1.GaussianPartChoiceSolve(EPS);
	A2.GaussianPartChoiceSolve(EPS);
	A3.GaussianPartChoiceSolve(EPS);
	A4.GaussianPartChoiceSolve(EPS);

	vector<vector<LT>> matrix_inv({ A1.SolutionX,A2.SolutionX, A3.SolutionX, A4.SolutionX });
	LT temp_el;
	for (int i = 0; i < n; ++i)
		for (int j = i + 1; j < n; ++j)
			if (i != j)
			{
				temp_el = matrix_inv[i][j];
				matrix_inv[i][j] = matrix_inv[j][i];
				matrix_inv[j][i] = temp_el;
			}
	return matrix_inv;
};