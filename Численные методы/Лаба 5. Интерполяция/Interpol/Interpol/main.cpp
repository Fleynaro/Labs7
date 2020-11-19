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

class LUDecomposition
{
	double** L, ** U;
	int* P;
	double sign;
public:
	int size;
	LUDecomposition(double** const matrix, const int size)
	{
		sign = 1.0;
		this->size = size;
		P = new int[size];
		L = new double* [size];
		U = new double* [size];
		for (size_t i = 0; i < size; i++)
		{
			L[i] = new double[size];
			U[i] = new double[size];
		}
		copyMatrixToMatrix(matrix, U, size);
		//иницилизируем подстановку P
		for (int i = 0; i < size; i++) {
			P[i] = i;
		}
		for (int k = 0; k < size; k++) {
			auto rowIdx = defineRowIdxWithMainValue(U, k, size);
			if (k != rowIdx) {
				//Смена строк
				swap(U[k], U[rowIdx]);
				swap(L[k], L[rowIdx]);
				swap(P[k], P[rowIdx]);
				sign *= -1.0;
			}

			//главный элемент
			double mainValue = U[k][k];
			//if (abs(mainValue) < minValue) {
			//	//Определяем ранг матрицы
			//	rank = k;
			//	return;
			//}
			//заполняем матрицу L
			for (int i = k; i < size; i++) {
				L[i][k] = U[i][k];
			}

			for (int j = k; j < size; j++) {
				U[k][j] /= mainValue;
			}
			//заполняем матрицу U
			for (int i = k + 1; i < size; i++) {
				for (int j = k; j < size; j++) {
					U[i][j] = U[i][j] - L[i][k] * U[k][j];
				}
			}
		}
	}
	~LUDecomposition()
	{
		for (size_t i = 0; i < size; i++)
		{
			delete[] L[i], U[i];
		}
		delete[] P, L, U;
	}
	//Решение уравнения Ax = b, то есть LUx = Pb
	void SolveSOLE(double* X, double* B) {
		//Ax = b (A = PLU)
		//LUx = Pb (Ux = y)
		//Ly = Pb
		double* vectorY = new double[size];
		//"умножаем" вектор b на матрицу перестановок
		double* vectorPB = new double[size];
		for (int i = 0; i < size; i++) {
			vectorPB[i] = B[P[i]];
		}
		SolveLy(L, vectorY, vectorPB, size);
		//printVector(Y, N);
		SolveUx(U, X, vectorY, size);
	}
	//Получение обратной матрицы
	void SolveBackwardMatrix(double** X) {
		//LUX = PE
		double* vectorX = new double[size];
		double* vectorE = new double[size];
		for (int t = 0; t < size; t++) {
			vectorE[t] = 0.0;
		}
		for (int i = 0; i < size; i++) {
			//формируем вектор Ei
			if (i != 0) {
				vectorE[i - 1] = 0.0;
			}
			vectorE[i] = 1.0;

			//получаем вектор-столбец X
			SolveSOLE(vectorX, vectorE);
			//записываем его в матрицу X
			for (int t = 0; t < size; t++) {
				X[t][i] = vectorX[t];
			}
		}
	}
private:
	//определить строку с главным элементом (который максимален в текущем столбце)
	int defineRowIdxWithMainValue(double** matrix, int k, int N) {
		int m = k;
		double maxValue = 0.0;
		for (int i = k; i < N; i++) {
			if (abs(matrix[i][k]) > maxValue) {
				m = i;
				maxValue = abs(matrix[m][k]);
			}
		}
		return m;
	}
	void copyMatrixToMatrix(double** srcMatrix, double** dstMatrix, int N) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				dstMatrix[i][j] = srcMatrix[i][j];
			}
		}
	}
	//Решение уравнения Ly = Pb
	void SolveLy(double** triangleMatrix, double* X, double* B, int N) {
		for (int i = 0; i < N; i++) {
			X[i] = B[i] / triangleMatrix[i][i];
			for (int j = 0; j < i; j++) {
				X[i] -= X[j] * triangleMatrix[i][j] / triangleMatrix[i][i];
			}
		}
	}
	//Решение уравнения Ux = y
	void SolveUx(double** triangleMatrix, double* X, double* B, int N) {
		for (int i = N - 1; i >= 0; i--) {
			X[i] = B[i];
			for (int j = N - 1; j > i; j--) {
				X[i] -= X[j] * triangleMatrix[i][j];
			}
		}
	}
};

class RmsApproximation
{
	using FuncType = std::function<double(double)>;
	double* m_fg; //столбец свободных членов
	double* m_c;
	double* m_f;
	double** m_gg; //матрица
	int m_vector_size;
	double m_step = 0.2;
	double normError = 0;
	double m_a, m_b; //концы отрезка
	FuncType m_gVector[3] = {
		[](double x) {
			return 1;
		},
		[](double x) {
			return x;
		},
		[](double x) {
			return x * x;
		}
	};
public:
	//Дискретный вариант
	RmsApproximation(FuncType function, double a, double b)
	{
		m_a = a;
		m_b = b;
		m_vector_size = 3;
		m_c = new double[m_vector_size];
		m_fg = new double[m_vector_size];
		//Свободные члены (f,g)
		for (size_t i = 0; i < m_vector_size; i++)
		{
			double sum = 0;
			for (double x = m_a; x <= m_b; x += m_step)
			{
				sum += function(x) * m_gVector[i](x);
			}
			m_fg[i] = sum;
		}
		//Создание матрицы
		m_gg = new double* [m_vector_size];
		for (size_t i = 0; i < m_vector_size; i++)
		{
			m_gg[i] = new double[m_vector_size];
		}
		//Ее заполнение (g,g)
		for (size_t i = 0; i < m_vector_size; i++)
		{
			for (size_t j = 0; j < m_vector_size; j++)
			{
				double sum = 0;
				for (double x = m_a; x <= m_b; x += m_step)
				{
					sum += m_gVector[i](x) * m_gVector[j](x);
				}
				m_gg[i][j] = sum;
			}
		}

		m_f = new double[1.0 / m_step + 1];
		double x = m_a;
		for (size_t i = 0; i < 1.0 / m_step + 1; i++)
		{
			m_f[i] = function(x);
			x += m_step;
		}
	}
	//Непрерывный вариант
	RmsApproximation(double gg[3][3], double fg[3])
	{
		m_a = 2;
		m_b = 1;
		m_vector_size = 3;
		m_c = new double[m_vector_size];
		m_fg = new double[m_vector_size];
		for (size_t i = 0; i < m_vector_size; i++)
		{
			m_fg[i] = fg[i];
		}
		m_gg = new double* [m_vector_size];
		for (size_t i = 0; i < m_vector_size; i++)
		{
			m_gg[i] = new double[m_vector_size];
		}
		for (size_t i = 0; i < m_vector_size; i++)
		{
			for (size_t j = 0; j < m_vector_size; j++)
			{
				m_gg[i][j] = gg[i][j];
			}
		}
	}
	~RmsApproximation()
	{
		for (size_t i = 0; i < m_vector_size; i++)
		{
			delete[] m_gg[i];
		}
		delete[] m_gg, m_fg, m_c, m_f;
	}
	void GetSolve()
	{
		LUDecomposition LUdec(m_gg, m_vector_size);
		LUdec.SolveSOLE(m_c, m_fg);
	}
	//Дискретный вариант
	void GetNormError()
	{
		//Погрешность
		FuncType fun_G = [&](double x) {
			return m_c[0] + (m_c[1] * x) + (m_c[2] * x * x);
		};
		double G = 0;
		for (double x = m_a; x < m_b; x += m_step)
		{
			G += pow(fun_G(x), 2);
		}
		double F = 0;
		for (size_t i = 0; i < 1.0 / m_step + 1; i++)
		{
			F += m_f[i] * m_f[i];
		}
		normError = sqrt(F - G);
	}
	//Непрерывный вариант
	void GetNormError(double ff, double gg)
	{
		normError = sqrt(abs(ff - gg));
	}
	void PrintSolve()
	{
		printf("Матрица\n");
		for (size_t i = 0; i < m_vector_size; i++)
		{
			for (size_t j = 0; j < m_vector_size; j++)
			{
				cout << fixed << setw(10) << setprecision(5) << m_gg[i][j] << "\t";
			}
			cout << endl;
		}
		printf("\nВектор правых частей\n");
		for (size_t i = 0; i < m_vector_size; i++)
		{
			printf("%.13f\t", m_fg[i]);
		}
		printf("\n\nP2(x) = (%.5f) + (%.5f)*x + (%.5f)*x^2\n", m_c[0], m_c[1], m_c[2]);
		printf("Норма погрешности: %.16f\n", normError);
	}

	double getIntegralOfG2(double x) {
		return pow(m_c[2], 2) * pow(x, 5) / 5 +
			m_c[1] * m_c[2] * pow(x, 4) / 2 +
			(pow(m_c[1], 2) + 2 * m_c[0] * m_c[2]) * pow(x, 3) / 3 +
			m_c[0] * m_c[1] * pow(x, 2) +
			pow(m_c[0], 2) * x;
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
		[](double x) { return pow(3, x) + 2 - x; },
		[](double x) { return pow(3, x) * log(3) - 1; }
	);

	//максимальная производная n-ого порядка функции f на отрезке [1, 2] (нужна для формулы оценки интерполяции)
	double M6 = 15.824; //6 порядка
	double M5 = 14.403; //5 порядка
	double M4 = 13.111; //4 порядка

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
		table.printTableRow({ x, func.m_func(x), Pn, abs(func.m_func(x) - Pn), abs(deltaEst) });
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
		double deltaEst2 = (M4 / 384 + M5 * h / 240) * pow(h, 4); //deltaEst < deltaEst2
		table3.printTableRow({ x, func.m_func(x), S31, abs(func.m_func(x) - S31), deltaEst });
	}

	printf("\nСреднеквадратичное приближение\n\n");
	RmsApproximation disRms(func.m_func, 1, 2);
	disRms.GetSolve();
	disRms.GetNormError();
	printf("Дискретный вариант\n");
	disRms.PrintSolve();

	double fg[3] = {
		5.96143535976102,
		9.349042367927889,
		15.148483258152980
	};
	double gg[3][3] = {
		{1., 3. / 2., 7. / 3.},
		{3. / 2., 7. / 3., 15. / 4.},
		{7. / 3., 15. / 4., 31. / 5.}
	};
	//double f[6] = { func.m_func(1.0), func.m_func(1.2), func.m_func(1.4), func.m_func(1.6), func.m_func(1.8), func.m_func(2.0), };
	RmsApproximation contRms(gg, fg);
	contRms.GetSolve();
	contRms.GetNormError(37.5829, contRms.getIntegralOfG2(B) - contRms.getIntegralOfG2(A));
	printf("\nНепрерывный вариант\n");
	contRms.PrintSolve();

	printf("\nРешение уравнения методом обратной интерполяции\n");
	//выбираем отрезок [a, b], где находится предполагаемый корень, а также выбираем n отрезков разбиения
	int n = 5;
	double a = 1.0;
	double b = 2.0;
	double c = func.m_func(1.5); //правая часть уравнения f(x) = c
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
