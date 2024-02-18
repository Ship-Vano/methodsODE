/*<<Чтение и запись из файла для конкретных приложений>>
 @author: Шаманов Иван
 2024*/
#pragma once
#include<fstream>
#include<vector>
#include<string>
using namespace std;
//TODO:чтение начальных условий из файла
/* Начальные условия задачи Коши
*	вектор u_0
*	[t_0, T]
*	f (id)
*	tau
*	tol
*/
/* Параметры метода
*	id метода
*/
//TODO:запись результатов в файл
template <typename DT>
void writeVectorToFile(ofstream& file, vector<DT> v)
{
	for (int i = 0; i < v.size(); ++i)
		file << v[i] << " ";
	file << " " << endl;
}
template<typename DT>
void write_data_to_file(string filepath, vector<vector<DT>> data)
{
	ofstream output_data;
	output_data.open(filepath);
	int n = data.size();
	int m = data[0].size();
	cout << m << endl;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			output_data << data[i][j] << "  ";
		}
		output_data << endl;
	}
	output_data.close();
}
/* Промежуточные решения
*	t0, y0
*	t1, y1
*	.....
*/