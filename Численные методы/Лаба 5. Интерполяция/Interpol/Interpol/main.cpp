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

//Разбиение отрезка [a, b]
struct Splitting
{
	int m_N;
	function<double(int)> m_nodes;
	Splitting(int N, function<double(int)> nodes)
		: m_N(N), m_nodes(nodes)
	{}

	//вычислить h(i)
	double getH(int i) {
		return m_nodes(i + 1) - m_nodes(i);
	}

	//вычислить индекс отрезка, которому принадлежит x
	int calcIndex(double x) {
		for (int i = 0; i < m_N; i++) {
			if (x >= m_nodes(i) && x <= m_nodes(i + 1))
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
	Splitting m_splitting;

	Interpolation(Function f, Splitting splitting)
		: m_f(f), m_splitting(splitting)
	{}

	//значение f(i)
	double getNodeY(int idx) {
		auto x = m_splitting.m_nodes(idx);
		return m_f.m_func(x);
	}

	//вычислить значение полученного полинома в произвольной точке
	virtual double interpolate(double x) = 0;
};

//Ньютоновская интерполяция
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
	double interpolate(double x) override {
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

//Кубическая сплайн интерполяция дефекта 1
class SplineInterpolation : public Interpolation
{
public:
	Vector m_vectorM;

	SplineInterpolation(Function f, Splitting splitting)
		: Interpolation(f, splitting)
	{
		m_vectorM = Vector(splitting.m_N + 1);
		m_vectorM[0] = f.m_derivative(m_splitting.m_nodes(0));
		m_vectorM[m_vectorM.m_N - 1] = f.m_derivative(m_splitting.m_nodes(m_vectorM.m_N - 1));
		calcVectorM();
	}

	//вычислить значение полученного полинома в произвольной точке
	double interpolate(double x) override {
		auto idx = m_splitting.calcIndex(x);
		auto h = m_splitting.getH(idx);
		auto tau = (x - m_splitting.m_nodes(idx)) / h;
		return Fi0(tau) * getNodeY(idx) + Fi0(1 - tau) * getNodeY(idx + 1) + h * (Fi1(tau) * m_vectorM[idx] - Fi1(1 - tau) * m_vectorM[idx + 1]);
	}

private:
	void calcVectorM() {
		Vector a(m_vectorM.m_N - 2);
		Vector c(a.m_N);
		Vector b(a.m_N);
		Vector f(a.m_N);
		Vector x(a.m_N);

		for (int i = 1, j = 0; i <= m_vectorM.m_N - 1; i++, j++) {
			auto curH = m_splitting.getH(i);
			auto prevH = m_splitting.getH(i - 1);

			b[j] = prevH / (prevH + curH);
			c[j] = 2;
			a[j] = 1 - b[j];
			f[j] = 3 * a[j] * (getNodeY(i) - getNodeY(i - 1)) / prevH +
				3 * b[j] * (getNodeY(i + 1) - getNodeY(i)) / curH;
		}

		f[0] -= m_splitting.getH(1) / (m_splitting.getH(0) + m_splitting.getH(1)) * m_vectorM[0];
		f[f.m_N - 1] -= m_splitting.getH(m_vectorM.m_N - 2) / (m_splitting.getH(m_vectorM.m_N - 2) + m_splitting.getH(m_vectorM.m_N - 1)) * m_vectorM[m_vectorM.m_N - 1];

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

	printf("\nИнтерполяция кубическим сплайном\n");
	SplineInterpolation splineInterpolation(func, splitting);

	printf("\nM5 = %f\n", M5);
	Table table2({ "x[i]", "df/dx(x[i])", "m[i]", "Delta", "Оценка" }, { 2, 7, 7, 15, 15 });
	table2.printTableHeader();
	for (int i = 0; i < N + 1; i++) {
		auto x = splitting.m_nodes(i);
		auto dF = func.m_derivative(x);
		auto m = splineInterpolation.m_vectorM[i];
		table2.printTableRow({ x, dF, m, abs(dF - m), abs(0.0) });
	}


	printf("\nM4 = %f\n", M4);
	Table table3({ "x", "f(x)", "S31(f;x)", "Abs(f(x)-S31(f;x))", "Оценка" }, { 2, 5, 5, 15, 15 });
	table3.printTableHeader();
	for (int i = 0; i < N; i++) {
		double x = 1.0 + (i + 0.5) * (2.0 - 1.0) / N;
		double S31 = splineInterpolation.interpolate(x);

		table3.printTableRow({ x, func.m_func(x), S31, abs(func.m_func(x) - S31), abs(0.0) });
	}


	/*Vector a(4);
	Vector b(4);
	Vector c(4);
	Vector x(4);
	Vector f(4);
	a[0] = 0;
	a[1] = 3;
	a[2] = 6;
	a[3] = 9;
	c[0] = 10;
	c[1] = 40;
	c[2] = 70;
	c[3] = 100;
	b[0] = 2;
	b[1] = 5;
	b[2] = 8;
	b[3] = 0;
	f[0] = 1;
	f[1] = 2;
	f[2] = 3;
	f[3] = 4;
	SplineInterpolation::SweepMethod(a, c, b, x, f);

	printf("\n\n");
	x.printVector();*/
}
