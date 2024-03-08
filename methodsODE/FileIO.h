/*<<������ � ������ �� ����� ��� ���������� ����������>>
 @author: ������� ����
 @author: ������ ����
 2024*/
#pragma once
#include<fstream>
#include<istream>
#include<vector>
#include<string>
#include<sstream>
using namespace std;
//TODO:������ ��������� ������� �� �����
/* ��������� ������� ������ ����
*	������ u_0
*	[t_0, T]
*	f (id)
*	tau
*	tol
*/
/* ��������� ������
*	id ������
*/
template<typename DT>
vector<vector<DT>> readInitVals(string path)
{
	/*string path = "InputData\\" + filename;*/
	vector<vector<DT>> result{};
	result.push_back({});
	string tmp_line;
	fstream input_data;
	input_data.open(path);
	if (input_data.is_open())
	{
		cout << "log[INFO]: Opening a file \""<< path << "\" to read..."<<endl;
		while(getline(input_data, tmp_line))
		{
			istringstream ss(tmp_line);
			result.emplace_back(istream_iterator<DT>(ss), \
				istream_iterator<double>());
		}
		result.erase(result.begin());
		input_data.close();
	}
	else
		cout << "log[ERROR]: Couldn't open or read the file "<< path << endl;
	return result;
}

//TODO:������ ����������� � ����
template <typename DT>
void writeVectorToFile(ofstream& file, vector<DT> v)
{
	for (int i = 0; i < v.size(); ++i)
		file << v[i] << " ";
	file << " " << endl;
}
template <typename DT>
void writeVectorToFile(fstream& file, vector<DT> v)
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
/* ������������� �������
*	t0, y0
*	t1, y1
*	.....
*/