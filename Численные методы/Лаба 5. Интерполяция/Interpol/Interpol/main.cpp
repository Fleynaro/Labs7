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

//Вектор
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
		for (int i = 0; i < N; i++) {
			m_vector[i] = 0;
		}
	}

	//обращение к компоненте вектора
	double& operator[](int index) {
		return m_vector[index];
	}

	//распечатка вектора
	void printVector(bool expFormat = false) {
		printf(" ");
		for (int i = 0; i < m_N; i++) {
			if (expFormat) {
				printf("%.10e ", m_vector[i]);
			}
			else {
				printf("%.7f ", m_vector[i]);
			}
		}
		printf("\n");
	}
};

//Матрица
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

//Функция и её производная
struct Function
{
	using FuncType = std::function<double(double)>;
	FuncType m_func;
	FuncType m_derivative;

	Function(FuncType func, FuncType derivative)
		: m_func(func), m_derivative(derivative)
	{}
};

//Таблица узлов интерполяции
class NodeTable
{
public:
	int m_N;
	function<double(int)> m_x_nodes;
	function<double(int)> m_y_nodes;

	NodeTable(int N, function<double(int)> x_nodes, function<double(int)> y_nodes)
		: m_N(N), m_x_nodes(x_nodes), m_y_nodes(y_nodes)
	{}

	//вычислить h(i)
	double getH(int i) {
		return m_x_nodes(i + 1) - m_x_nodes(i);
	}

	//вычислить индекс отрезка, которому принадлежит x
	int calcIndex(double x) {
		for (int i = 0; i < m_N; i++) {
			if (x >= m_x_nodes(i) && x <= m_x_nodes(i + 1))
				return i;
		}
		throw std::logic_error("x out of the range");
	}
};

//Интерполяция (абстрактный класс)
class Interpolation
{
public:
	Function m_f;
	NodeTable m_nodeTable;

	Interpolation(Function f, NodeTable splittingTable)
		: m_f(f), m_nodeTable(splittingTable)
	{}

	//вычислить значение полученного полинома в произвольной точке
	virtual double interpolate(double x) = 0;
};

//Ньютоновская интерполяция
class NewtonianInterpolation : public Interpolation
{
	//таблица разделенных разностей (без столбца значений X)
	Matrix m_divDifferTable;
public:
	NewtonianInterpolation(Function f, NodeTable nodeTable)
		: Interpolation(f, nodeTable)
	{
		m_divDifferTable = Matrix(m_nodeTable.m_N + 1);
		fillDivDifferTable();
	}

	//вычислить значение полученного полинома в произвольной точке
	double interpolate(double x) override {
		double w = 1;
		double result = m_divDifferTable[0][0];
		for (int i = 1; i < m_divDifferTable.m_N; i++) {
			w *= x - m_nodeTable.m_x_nodes(i - 1);
			result += m_divDifferTable[0][i] * w;
		}
		return result;
	}

	//распечатка таблицы разделенных разностей
	void printDivDifferTable()
	{
		printf("Таблица разделенных разностей\n");
		for (int i = 0; i < m_divDifferTable.m_N; i++) {
			printf("%.2f  ", m_nodeTable.m_x_nodes(i));
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
			m_divDifferTable[i][0] = m_nodeTable.m_y_nodes(i);
		}

		//заполняем все последующие столбцы
		for (int j = 1; j < m_divDifferTable.m_N; j++) {
			for (int i = 0; i < m_divDifferTable.m_N - j; i++) {
				m_divDifferTable[i][j] =
					(m_divDifferTable[i + 1][j - 1] - m_divDifferTable[i][j - 1]) / (m_nodeTable.m_x_nodes(i + j) - m_nodeTable.m_x_nodes(i));
			}
		}
	}
};

//Кубическая сплайн интерполяция дефекта 1
class SplineInterpolation : public Interpolation
{
public:
	Vector m_vectorM;

	SplineInterpolation(Function f, NodeTable nodeTable)
		: Interpolation(f, nodeTable)
	{
		m_vectorM = Vector(m_nodeTable.m_N + 1);
		m_vectorM[0] = f.m_derivative(m_nodeTable.m_x_nodes(0));
		m_vectorM[m_vectorM.m_N - 1] = f.m_derivative(m_nodeTable.m_x_nodes(m_vectorM.m_N - 1));
		calcVectorM();
	}

	//вычислить значение полученного полинома в произвольной точке
	double interpolate(double x) override {
		return interpolate(x, nullptr, 0, 0);
	}

	//интерполяция + вычисление погрешности
	double interpolate(double x, double* errorEst, double M4, double M5) {
		auto idx = m_nodeTable.calcIndex(x);
		auto h = m_nodeTable.getH(idx);
		auto tau = (x - m_nodeTable.m_x_nodes(idx)) / h;
		if (errorEst)
			*errorEst = tau * Fi1(tau) * pow(h, 4) * M4 / 24 + (Fi1(tau) + Fi1(1 - tau)) * pow(h, 5) * M5 / 60;
		return Fi0(tau) * m_nodeTable.m_y_nodes(idx) + Fi0(1 - tau) * m_nodeTable.m_y_nodes(idx + 1) + h * (Fi1(tau) * m_vectorM[idx] - Fi1(1 - tau) * m_vectorM[idx + 1]);
	}

private:
	//найти значения m(i)
	void calcVectorM() {
		Vector a(m_vectorM.m_N - 2);
		Vector c(a.m_N);
		Vector b(a.m_N);
		Vector f(a.m_N);
		Vector x(a.m_N);

		for (int i = 1, j = 0; i <= m_vectorM.m_N - 1; i++, j++) {
			auto curH = m_nodeTable.getH(i);
			auto prevH = m_nodeTable.getH(i - 1);

			b[j] = prevH / (prevH + curH);
			c[j] = 2;
			a[j] = 1 - b[j];
			f[j] = 3 * a[j] * (m_nodeTable.m_y_nodes(i) - m_nodeTable.m_y_nodes(i - 1)) / prevH +
				3 * b[j] * (m_nodeTable.m_y_nodes(i + 1) - m_nodeTable.m_y_nodes(i)) / curH;
		}

		f[0] -= m_nodeTable.getH(1) / (m_nodeTable.getH(0) + m_nodeTable.getH(1)) * m_vectorM[0];
		f[f.m_N - 1] -= m_nodeTable.getH(m_vectorM.m_N - 2) / (m_nodeTable.getH(m_vectorM.m_N - 2) + m_nodeTable.getH(m_vectorM.m_N - 1)) * m_vectorM[m_vectorM.m_N - 1];

		SweepMethod(a, c, b, x, f);
		for (int i = 0; i < x.m_N; i++) {
			m_vectorM[i + 1] = x[i];
		}
	}

	double Fi0(double tau) {
		return (1 + 2 * tau) * pow(1 - tau, 2);
	}

	double Fi1(double tau) {
		return tau * pow(1 - tau, 2);
	}

	//Метод прогонки (матрица состоит из 3-х диагоналей: a, c, b)
	static void SweepMethod(Vector& a, Vector& c, Vector& b, Vector& x, Vector& f) {
		Vector alpha(x.m_N);
		Vector beta(x.m_N);
		alpha[0] = -b[0] / c[0]; //b1 / c1
		beta[0] = f[0] / c[0]; //f1 / c1

		for (int i = 0, j = 1; i < alpha.m_N - 1; i++, j++) {
			double denominator = c[j] + a[j] * alpha[i]; //c(j) - a(j) * alpha(j)
			alpha[i + 1] = -b[j] / denominator; //b(j) / den
			beta[i + 1] = (f[j] - a[j] * beta[i]) / denominator; //[f(j) + a(j) * beta(j)] / den
		}

		x[x.m_N - 1] = beta[beta.m_N - 1];
		for (int i = x.m_N - 2; i >= 0; i--) {
			x[i] = alpha[i] * x[i + 1] + beta[i];
		}
	}
};

int main()
{
	system("chcp 1251");

	int N = 5;
	double A = 1;
	double B = 2;

	//функция и её производная
	Function func(
		[](double x) { return pow(3, x - 1) + 4 - x; },
		[](double x) { return pow(3, x - 1) * log(3) - 1; }
	);

	//максимальная производная n-ого порядка функции f на отрезке [1, 2] (нужна для формулы оценки интерполяции)
	double M6 = 5.2745810; //6 порядка
	double M5 = 4.8011306; //5 порядка
	double M4 = 4.3701774; //4 порядка

	printf("Интерполяционная формула Ньютона\n");
	auto nodesX = [&](int i) { return A + i * (B - A) / N; };
	auto nodesY = [&](int i) { return func.m_func(nodesX(i)); };
	NodeTable nodeTable(N, nodesX, nodesY);
	NewtonianInterpolation newtonianInterpolation(func, nodeTable);
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
			w *= x - nodeTable.m_x_nodes(i);
		double deltaEst = w * M6 / 720;
		table.printTableRow({ x, func.m_func(x), Pn, abs(func.m_func(x)- Pn), abs(deltaEst) });
	}
	
	printf("\nИнтерполяция кубическим сплайном\n");
	SplineInterpolation splineInterpolation(func, nodeTable);

	//Вычисление производных в узлах интерполяции
	printf("\nM5 = %f\n", M5);
	Table table2({ "x[i]", "df/dx(x[i])", "m[i]", "Delta", "Оценка" }, { 2, 7, 7, 15, 15 });
	table2.printTableHeader();
	for (int i = 0; i < N + 1; i++) {
		auto x = nodeTable.m_x_nodes(i);
		auto dF = func.m_derivative(x);
		auto m = splineInterpolation.m_vectorM[i];
		double h = nodeTable.getH(0);
		double deltaEst = M5 * pow(h, 4) / 60;
		table2.printTableRow({ x, dF, m, abs(dF - m), deltaEst });
	}

	//Вычисление промежуточных значений
	printf("\nM4 = %f\n", M4);
	Table table3({ "x", "f(x)", "S31(f;x)", "Abs(f(x)-S31(f;x))", "Оценка" }, { 2, 5, 5, 15, 15 });
	table3.printTableHeader();
	for (int i = 0; i < N; i++) {
		double x = 1.0 + (i + 0.5) * (2.0 - 1.0) / N;
		double deltaEst;
		double S31 = splineInterpolation.interpolate(x, &deltaEst, M4, M5);
		double h = nodeTable.getH(0);
		double deltaEst2 = (M4 / 384 + M5 * h / 240)* pow(h, 4); //deltaEst < deltaEst2
		table3.printTableRow({ x, func.m_func(x), S31, abs(func.m_func(x) - S31), deltaEst });
	}




	//2 пункт лабы ...




	printf("\nРешение уравнения методом обратной интерполяции\n");
	//выбираем отрезок [a, b], где находится предполагаемый корень, а также выбираем n отрезков разбиения
	int n = 5;
	double a = 1.0;
	double b = 2.0;
	double c = 4.2321; //правая часть уравнения f(x) = c
	auto nodesX2 = [&](int i) { return a + i * (b - a) / n; };
	auto nodesY2 = [&](int i) { return func.m_func(nodesX2(i)); };
	NodeTable nodeTable2(n, nodesY2, nodesX2); //наоборот
	NewtonianInterpolation newtonianInterpolation2(func, nodeTable2);
	newtonianInterpolation2.printDivDifferTable();
	auto x = newtonianInterpolation2.interpolate(c);

	printf("c = %.4f\n", c);
	printf("Корень=%.5f\n", x);
	printf("Невязка=Abs(f(x)-c)=%f\n", abs(func.m_func(x) - c));
}
