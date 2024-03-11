/*<<Методы решений Нелинейных уравнений>>
 @author: Шаманов Иван
 @author: Климов Олег
 2023*/

#pragma once
#include<iostream>
#include<map>
#include"LinOp.h"
#include"GaussianSolve.h"

using namespace std;

template<typename DT, typename F>
map<DT, DT> UniformNet(F func, DT a, DT b, int n)
{
	map<DT, DT> net;
	vector<DT> xs;
	vector<DT> ys;
	DT h = fabs((b - a) / (n - 1));
	for (size_t i = 0; i < n; ++i)
	{
		xs.push_back(a + i * h);
		ys.push_back(func(xs[i]));
		net[xs[i]] = ys[i];
	}
	return net;
}

template<typename DT, typename F>
vector<vector<DT>> RootLocalization(F func, DT a, DT b, int N)
{
	int rootAmnt = 0;
	vector<vector<DT>> rootAreas{};
	map<DT, DT> table = UniformNet(func, a, b, N);
	//Проходимся циклом по всем парам, находим число корней:
	for (auto it = std::next(table.begin()); it != table.end(); ++it) {
		auto prev_el = std::prev(it);
		DT y_now = it->second;
		DT y_prev = prev_el->second;
		if (y_now * y_prev < 0)
		{
			++rootAmnt;
			rootAreas.push_back({ prev_el->first, it->first });
		}
	}
	//cout << rootAmnt<<endl;
	//out(rootAreas);
	//Измельчаем сетку
	//Проходимся по всем парам, находим число корней
	//Если одинаковые числа корней, то оставляем
	//Иначе: производим ещё измельчение и последующую проверку и тд...
	return rootAreas;
}

template<typename DT, typename F>
DT BisectionSolve(F func, DT a, DT b, DT eps = 1e-6, int max_its = 1000)
{
	DT f1 = func(a);
	DT f2 = func(b);
	if (f1 * f2 > 0) return NULL;
	/*if (abs(f1) < eps)
	{
		return a;
	}
	if (abs(f2) < eps)
	{
		return b;
	}*/
	//для определённости f1<=0, а f2 >=0
	int sign = 1;
	if (f1 > 0 && f2 < 0) sign = -1;
	int its = 0;
	DT c_prev = 0;
	while (its < max_its)
	{
		DT c = (a + b) / 2;
		if (abs(c - c_prev) < eps)
		{
			cout << "Число итераций: " << its << endl;
			return c;
		}
		DT f_mid = sign * func(c);
		/*if (abs(f_mid)<eps)
		{
			cout << "Число итераций: " << its << endl;
			return (a + b) / 2;
		}*/
		if (f_mid > 0)
		{
			b = c;
		}
		if (f_mid < 0)
		{
			a = c;
		}
		if ((b - a) / 2 < eps)
		{
			cout << "Число итераций: " << its << endl;
			return c;
		}
		c_prev = c;
		++its;
	}
	cout << "Число итераций: " << its << endl;
	return NULL;
}

//М. Ньютона с численным вычислением производной
template<typename DT, typename F>
DT deriv(F func, DT x, DT eps = 1e-6)
{
	DT h = eps;
	return (func(x + h) - func(x - h)) / (2 * h);
}
template<typename DT, typename F>
DT NewtonSolve(F func, DT a, DT b, DT eps = 1e-6, int max_its = 1000, DT x_0 = (a + b) / 2)
{
	cout << x_0 << endl;
	DT x_k = x_0;
	int its = 0;
	vector<DT> xs{};
	while (its < max_its)
	{
		xs.push_back(x_k);
		DT x_prev = x_k;
		DT f_k = func(x_k);
		/*if (abs(f_k) < eps)
		{
			cout << "Число итераций: " << its << endl;
			return x_k;
		}*/
		DT df_k = deriv(func, x_k);
		x_k = x_prev - f_k / df_k;
		if (abs(x_k - x_prev) < eps)
		{
			cout << "Число итераций: " << its << endl;
			return x_k;
		}
		///////проверка на выход из области и устранение
		if ((x_k < a) || (x_k > b))
		{
			DT reduce_coef = 0.5;
			DT alpha = 0.5;
			int inner_its = 1;
			x_k = x_prev - alpha * f_k / df_k;
			while (((x_k < a) && (x_k > b)) && (inner_its < max_its))
			{
				alpha = alpha * reduce_coef;
				x_k = x_prev - alpha * f_k / df_k;
				++inner_its;
			}
		}
		//////
		//////Проверка на "зацикливание" и устранение
		for (auto x_known : xs)
		{
			if (abs(x_k - x_known) < eps * 1e-6)
			{
				//Если новое приближение почти совпадает с одним из
				//ранее найденных, то вносим небольшую погрешность
				if (abs(b - x_k) < eps)
					x_k = x_k - eps;
				else
					x_k = x_k + eps;
				xs.clear();
				break;
			}
		}
		//////
		++its;
	}
	cout << "Число итераций: " << its << endl;
	return NULL;
}

//М. Ньютона с аналитическим вычислением производной
template<typename DT, typename F, typename dF>
DT NewtonSolve(F func, dF dfunc, DT a, DT b, DT eps = 1e-6, int max_its = 1000, DT x_0 = (a + b) / 2)
{
	DT x_k = x_0;
	int its = 0;
	vector<DT> xs{};
	while (its < max_its)
	{
		xs.push_back(x_k);
		DT x_prev = x_k;
		DT f_k = func(x_k);
		/*if (abs(f_k) < eps)
		{
			cout << "Число итераций: " << its << endl;
			return x_k;
		}*/
		DT df_k = dfunc(x_k);
		x_k = x_prev - f_k / df_k;
		if (abs(x_k - x_prev) < eps)
		{
			cout << "Число итераций: " << its << endl;
			return x_k;
		}
		///////проверка на выход из области и устранение
		if ((x_k < a) || (x_k > b))
		{
			DT reduce_coef = 0.5;
			DT alpha = 0.5;
			int inner_its = 1;
			x_k = x_prev - alpha * f_k / df_k;
			while (((x_k < a) && (x_k > b)) && (inner_its < max_its))
			{
				alpha = alpha * reduce_coef;
				x_k = x_prev - alpha * f_k / df_k;
				++inner_its;
			}
		}
		//////
		//////Проверка на "зацикливание" и устранение
		for (auto x_known : xs)
		{
			if (abs(x_k - x_known) < eps * 1e-6)
			{
				//Если новое приближение почти совпадает с одним из
				//ранее найденных, то вносим небольшую погрешность
				if (abs(b - x_k) < eps)
					x_k = x_k - eps;
				else
					x_k = x_k + eps;
				xs.clear();
				break;
			}
		}
		//////
		++its;
	}
	cout << "Число итераций: " << its << endl;
	return NULL;
}

template<typename DT, typename F>
vector<DT> Solve(F func, DT a, DT b, int N = 1000, DT eps = 0.01, int max_its = 10000, string method = "Bisection", DT x_0 = NULL)
{
	vector<vector<DT>> rootAreas = RootLocalization(func, a, b, N);
	vector<DT> roots{};
	if (method == "Bisection")
		for (int i = 0; i < rootAreas.size(); ++i)
		{
			double leftBorder = rootAreas[i][0];
			double rightBorder = rootAreas[i][1];
			//cout << "---l: " << leftBorder << "  r: " << rightBorder << endl;
			DT root = BisectionSolve(func, leftBorder, rightBorder, eps);
			if (root != NULL)
				roots.push_back(root);
		}
	if (method == "Newton")
		for (int i = 0; i < rootAreas.size(); ++i)
		{
			double leftBorder = rootAreas[i][0];
			double rightBorder = rootAreas[i][1];
			//cout << "---l: " << leftBorder << "  r: " << rightBorder << endl;
			DT root = NewtonSolve(func, leftBorder, rightBorder, eps, max_its, x_0);
			if (root != NULL)
				roots.push_back(root);
		}
	return roots;
}
template<typename DT, typename F, typename dF>
vector<DT> Solve(F func, dF dfunc, DT a, DT b, int N = 1000, DT eps = 0.01, int max_its = 10000, string method = "Bisection", DT x_0 = NULL)
{
	vector<vector<DT>> rootAreas = RootLocalization(func, a, b, N);
	vector<DT> roots{};
	if (method == "Bisection")
		for (int i = 0; i < rootAreas.size(); ++i)
		{
			double leftBorder = rootAreas[i][0];
			double rightBorder = rootAreas[i][1];
			//cout << "---l: " << leftBorder << "  r: " << rightBorder << endl;
			DT root = BisectionSolve(func, leftBorder, rightBorder, eps);
			if (root != NULL)
				roots.push_back(root);
		}
	if (method == "Newton")
		for (int i = 0; i < rootAreas.size(); ++i)
		{
			double leftBorder = rootAreas[i][0];
			double rightBorder = rootAreas[i][1];
			//cout << "---l: " << leftBorder << "  r: " << rightBorder << endl;
			DT root = NewtonSolve(func, dfunc, leftBorder, rightBorder, eps, max_its, x_0);
			if (root != NULL)
				roots.push_back(root);
		}
	return roots;
}
//М. Ньютона для систем с численным вычислением производных
template<typename DT, typename VT, typename F>
DT deriv_x(F func, VT x, DT eps = 1e-6)
{
	DT h = eps;
	return (func({ x[0] + h, x[1] }) - func({ x[0], x[1] })) / h;
	//return (func({ x[0] + h, x[1]}) - func({x[0] - h, x[1]})) / (2 * h);
}
template<typename DT, typename VT, typename F>
DT deriv_y(F func, VT x, DT eps = 1e-6)
{
	DT h = eps;
	return (func({ x[0], x[1] + h }) - func({ x[0], x[1] })) / h;
	//return (func({ x[0], x[1] + h }) - func({ x[0], x[1] - h })) / (2 * h);
}

template<typename DT, typename F1, typename F2>
vector<vector<DT>> JacobyMatrix(F1 func1, F2 func2, vector<DT> x, DT eps = 1e-6)
{
	DT a = deriv_x(func1, x, eps);
	DT b = deriv_y(func1, x, eps);
	DT c = deriv_x(func2, x, eps);
	DT d = deriv_y(func2, x, eps);
	vector<vector<DT>> result{ { a, b }, {c, d} };
	return result;
}

template<typename DT, typename F1, typename F2>
vector<DT> NewtonSolve(F1 func1, F2 func2, vector<DT> range1, vector<DT> range2, DT eps = 1e-6, int max_its = 1000, vector<DT> x_0 = NULL, bool LOG = false)
{
	int its = 0;
	vector<DT> x_k = x_0;
	while (its < max_its)
	{
		vector<DT> x_prev = x_k;
		vector<DT> f_k{ func1(x_k), func2(x_k) };
		vector<vector<DT>> J_k = JacobyMatrix(func1, func2, x_k, eps);
		DT det_J = J_k[0][0] * J_k[1][1] - J_k[0][1] * J_k[1][0];
		vector<vector<DT>> J_k_inv{ {J_k[1][1], -J_k[0][1]},
								   {-J_k[1][0], J_k[0][0]} };
		J_k_inv = (1 / det_J) * J_k_inv;
		vector<DT> product = J_k_inv * f_k;
		x_k = x_prev + (-1) * product;
		++its;
		///////проверка на выход из области и устранение
		if ((x_k[0] < range1[0]) || (x_k[0] > range1[1]) || \
			(x_k[1] < range2[0]) || (x_k[1] > range2[1]))
		{
			DT reduce_coef = 0.5;
			DT alpha = 0.5;
			int inner_its = 1;
			x_k = x_prev + (-1) * alpha * product;
			while (((x_k[0] < range1[0]) || (x_k[0] > range1[1]) || \
				(x_k[1] < range2[0]) || (x_k[1] > range2[1])) && (inner_its < max_its))
			{
				alpha = alpha * reduce_coef;
				x_k = x_prev + (-1) * alpha * product;
				++inner_its;
			}
		}
		//////
		if (vec_norm_inf(x_k - x_prev) < eps)
		{
			break;
		}
	}
	if (LOG) x_k.push_back(its);
	return x_k;
}

template<typename DT, typename VT, typename F>
DT deriv(F func, VT x, int arg_pos, int f_pos=0, DT eps = 1e-6)
{
	DT h = eps;
	VT xh = x;
	xh[arg_pos] += h;
	vector<DT> f_0 = func(x);
	vector<DT> f_1 = func(xh);
	int dimension_size = f_0.size();
	if (dimension_size > 1)
	{
		return (1/h)*(f_1[f_pos] - f_0[f_pos]) ;
	}
}
template<typename DT, typename F>
vector<vector<DT>> JacobyMatrix(F func, vector<DT> x, int dim_size = 2,DT eps = 1e-6)
{
	/*DT a = deriv_x(func1, x, eps);
	DT b = deriv_y(func1, x, eps);
	DT c = deriv_x(func2, x, eps);
	DT d = deriv_y(func2, x, eps);
	vector<vector<DT>> result{ { a, b }, {c, d} };
	return result;*/
	vector<vector<DT>> result;
	if (dim_size == 2)
	{
		DT a = deriv(func, x, 0, 0, eps);
		DT b = deriv(func, x, 1, 0, eps);
		DT c = deriv(func, x, 0, 1, eps);
		DT d = deriv(func, x, 1, 1, eps);
		result = { { a, b }, { c, d } };
	}
	else if (dim_size == 3)
	{
		/*
		  x0 x1 x2
		f0 a b c 
		f1 d e f
	    f2 g h i
		*/
		DT a = deriv(func, x, 0, 0, eps);
		DT b = deriv(func, x, 1, 0, eps);
		DT c = deriv(func, x, 2, 0, eps);
		DT d = deriv(func, x, 0, 1, eps);
		DT e = deriv(func, x, 1, 1, eps);
		DT f = deriv(func, x, 2, 1, eps);
		DT g = deriv(func, x, 0, 2, eps);
		DT h = deriv(func, x, 1, 2, eps);
		DT i = deriv(func, x, 2, 2, eps);
		result = { {a, b, c}, {d, e, f}, {g, h, i} };
	}
	return result;
}
template<typename DT, typename F>
vector<DT> NewtonSolve(F func, vector<DT> range, DT eps = 1e-6, int max_its = 1000, vector<DT> x_0 = NULL, bool LOG = false)
{
	int its = 0;
	vector<DT> x_k = x_0;
	while (its < max_its)
	{
		vector<DT> x_prev = x_k;
		vector<DT> f_k = func(x_k);
		int dim_size = f_k.size();
		vector<vector<DT>> J_k = JacobyMatrix(func, x_k, dim_size, eps);
		vector<vector<DT>> J_k_inv = InverseMatrix(J_k, eps);
		vector<DT> product = J_k_inv * f_k;
		x_k = x_prev + (-1) * product;
		++its;
		///////проверка на выход из области и устранение
		if ((x_k[0] < range[0]) || (x_k[0] > range[1]) || \
			(x_k[1] < range[0]) || (x_k[1] > range[1]))
		{
			DT reduce_coef = 0.5;
			DT alpha = 0.5;
			int inner_its = 1;
			x_k = x_prev + (-1) * alpha * product;
			while (((x_k[0] < range[0]) || (x_k[0] > range[1]) || \
				(x_k[1] < range[0]) || (x_k[1] > range[1])) && (inner_its < max_its))
			{
				alpha = alpha * reduce_coef;
				x_k = x_prev + (-1) * alpha * product;
				++inner_its;
			}
		}
		//////
		if (vec_norm_inf(x_k - x_prev) < eps)
		{
			break;
		}
	}
	if (LOG) x_k.push_back(its);
	return x_k;
}

//М. Ньютона для систем с аналитическим вычислением производных
template<typename DT, typename dF1x, typename dF1y, \
	typename dF2x, typename dF2y>
	vector<vector<DT>> JacobyMatrix(dF1x df1x, dF1y df1y, \
		dF2x df2x, dF2y df2y, vector<DT> x)
{
	DT a = df1x(x);
	DT b = df1y(x);
	DT c = df2x(x);
	DT d = df2y(x);
	vector<vector<DT>> result{ { a, b }, {c, d} };
	return result;
}
template<typename DT, typename F1, typename F2, \
	typename dF1x, typename dF1y, typename dF2x, typename dF2y>
	vector<DT> NewtonSolve(F1 func1, F2 func2, dF1x df1x, dF1y df1y, \
		dF2x df2x, dF2y df2y, vector<DT> range1, vector<DT> range2, DT eps = 1e-6, int max_its = 1000, vector<DT> x_0 = NULL, bool LOG = false)
{
	int its = 0;
	vector<DT> x_k = x_0;
	while (its < max_its)
	{
		vector<DT> x_prev = x_k;
		vector<DT> f_k{ func1(x_k), func2(x_k) };
		vector<vector<DT>> J_k = JacobyMatrix(df1x, df1y, df2x, df2y, x_k);
		DT det_J = J_k[0][0] * J_k[1][1] - J_k[0][1] * J_k[1][0];
		vector<vector<DT>> J_k_inv{ {J_k[1][1], -J_k[0][1]},
								   {-J_k[1][0], J_k[0][0]} };
		J_k_inv = (1 / det_J) * J_k_inv;
		vector<DT> product = J_k_inv * f_k;
		x_k = x_prev + (-1) * product;
		++its;
		///////проверка на выход из области и устранение
		if ((x_k[0] < range1[0]) || (x_k[0] > range1[1]) || \
			(x_k[1] < range2[0]) || (x_k[1] > range2[1]))
		{
			DT reduce_coef = 0.5;
			DT alpha = 0.5;
			int inner_its = 1;
			x_k = x_prev + (-1) * alpha * product;
			while (((x_k[0] < range1[0]) || (x_k[0] > range1[1]) || \
				(x_k[1] < range2[0]) || (x_k[1] > range2[1])) && (inner_its < max_its))
			{
				alpha = alpha * reduce_coef;
				x_k = x_prev + (-1) * alpha * product;
				++inner_its;
			}
		}
		//////
		if (vec_norm_inf(x_k - x_prev) < eps)
			break;
	}
	//cout << "Число итераций: " << its << endl;
	if (LOG) x_k.push_back(its);
	return x_k;
}