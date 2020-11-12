#include <iostream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <functional>
#include <iomanip>
#include <cstdarg>

using namespace std;

//Таблица
class Table
{
	list<string> m_colnames;
	list<int> m_precisions;
public:
	Table(const list<string>& colnames, const list<int>& precisions)
		: m_colnames(colnames), m_precisions(precisions)
	{}

	//распечатка заголовка таблицы
	void printTableHeader()
	{
		cout << "|";
		auto it = m_precisions.begin();
		for (const auto& colname : m_colnames)
		{
			cout << setw(*it * 2) << colname << "|";
			it++;
		}
		cout << "\n";
	}

	//распечатка строки таблицы
	void printTableRow(const vector<double>& values)
	{
		cout << "|";
		auto it = m_precisions.begin();
		for (const auto& value : values)
		{
			cout << setw(*it * 2) << setprecision(*it) << fixed << value << "|";
			it++;
		}
		cout << "\n";
	}
};

class Vector
{
	double* m_vector;
public:
	int m_N;

	Vector() = default;

	Vector(int N)
		: m_N(N)
	{
		m_vector = new double[N];
	}

	double& operator[](int index) {
		return m_vector[index];
	}
};

class Matrix
{
	Vector* m_vectors;
public:
	int m_N;

	Matrix() = default;

	Matrix(int N)
		: m_N(N)
	{
		m_vectors = new Vector[N];
		for (int i = 0; i < N; i++) {
			m_vectors[i] = Vector(N);
		}
	}

	Vector& operator[](int index) {
		return m_vectors[index];
	}
};

//Функция
struct Function
{
	using FuncType = std::function<double(double)>;
	FuncType m_func;
	FuncType m_derivative;

	Function(FuncType func, FuncType derivative)
		: m_func(func), m_derivative(derivative)
	{}
};

struct Splitting
{
	int m_N;
	function<double(int)> m_nodes;
	Splitting(int N, function<double(int)> nodes)
		: m_N(N), m_nodes(nodes)
	{}
};

class Interpolation
{
public:
	Function m_f;
	Splitting m_splitting;

	Interpolation(Function f, Splitting splitting)
		: m_f(f), m_splitting(splitting)
	{}

	double getNodeY(int idx) {
		auto x = m_splitting.m_nodes(idx);
		return m_f.m_func(x);
	}
};

class NewtonianInterpolation : public Interpolation
{
	//таблица разделенных разностей (без столбца значений X)
	Matrix m_divDifferTable;
public:
	NewtonianInterpolation(Function f, Splitting splitting)
		: Interpolation(f, splitting)
	{
		m_divDifferTable = Matrix(splitting.m_N + 1);
		fillDivDifferTable();
	}

	//вычислить значение полученного полинома в произвольной точке
	double interpolate(double x) {
		double w = 1;
		double result = m_divDifferTable[0][0];
		for (int i = 1; i < m_divDifferTable.m_N; i++) {
			w *= x - m_splitting.m_nodes(i - 1);
			result += m_divDifferTable[0][i] * w;
		}
		return result;
	}

	//распечатка таблицы разделенных разностей
	void printDivDifferTable()
	{
		printf("Таблица разделенных разностей\n");
		for (int i = 0; i < m_divDifferTable.m_N; i++) {
			printf("%.2f  ", m_splitting.m_nodes(i));
			for (int j = 0; j < m_divDifferTable.m_N - i; j++) {
				printf("%.6f  ", m_divDifferTable[i][j]);
			}
			printf("\n");
		}
	}
private:
	void fillDivDifferTable() {
		//заполняем первый столбец
		for (int i = 0; i < m_divDifferTable.m_N; i++) {
			m_divDifferTable[i][0] = getNodeY(i);
		}

		//заполняем все последующие столбцы
		for (int j = 1; j < m_divDifferTable.m_N; j++) {
			for (int i = 0; i < m_divDifferTable.m_N - j; i++) {
				m_divDifferTable[i][j] =
					(m_divDifferTable[i + 1][j - 1] - m_divDifferTable[i][j - 1]) / (m_splitting.m_nodes(i + j) - m_splitting.m_nodes(i));
			}
		}
	}
};

class SplineInterpolation : public Interpolation
{

};

int main()
{
	system("chcp 1251");

	int N = 5;

	//функция и её производная
	Function func(
		[](double x) { return pow(3, x - 1) + 4 - x; },
		[](double x) { return pow(3, x - 1) * log(3) - 1; }
	);

	//максимальная производная (нужна для формулы оценки интерполяции)
	double M6 = 5.2745810; //для ньютона
	double M5 = 4.8011306; //для сплайнов
	double M4 = 4.3701774;

	printf("Интерполяционная формула Ньютона\n");
	Splitting splitting(N, [N](int i) { return 1.0 + i * (2.0 - 1.0) / N; });
	NewtonianInterpolation newtonianInterpolation(func, splitting);
	newtonianInterpolation.printDivDifferTable();

	printf("\nM6 = %f\n", M6);
	Table table({ "x", "f(x)", "Pn(x)", "Delta", "Оценка" }, { 2, 5, 5, 15, 15 });
	table.printTableHeader();

	for (int i = 0; i < N; i++) {
		double x = 1.0 + (i + 0.5) * (2.0 - 1.0) / N;
		double Pn = newtonianInterpolation.interpolate(x);

		//расчет оценки
		double w = 1;
		for (int i = 0; i <= N; i++)
			w *= x - splitting.m_nodes(i);
		double deltaEst = w * M6 / 720;
		table.printTableRow({ x, func.m_func(x), Pn, abs(func.m_func(x)- Pn), abs(deltaEst) });
	}



}
