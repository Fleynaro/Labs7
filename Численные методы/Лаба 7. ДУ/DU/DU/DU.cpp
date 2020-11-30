#include <iostream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <functional>
#include <iomanip>
#include <cstdarg>
#include <list>
#include <vector>

using namespace std;

//Таблица
class ConsoleTable
{
	list<string> m_colnames;
	list<int> m_precisions;
	bool m_itr_col;
public:
	ConsoleTable(const list<string>& colnames, const list<int>& precisions, bool itr_col)
		: m_colnames(colnames), m_precisions(precisions), m_itr_col(itr_col)
	{}

	//распечатка заголовка таблицы
	void printTableHeader()
	{
		cout << "|";
		if (m_itr_col) {
			cout << setw(4) << "Itr |";
		}
		auto it = m_precisions.begin();
		for (const auto& colname : m_colnames)
		{
			cout << setw(*it * 2) << colname << "|";
			it++;
		}
		cout << "\n";
	}

	//распечатка строки таблицы
	void printTableRow(int N, const vector<double>& values)
	{
		cout << "|";
		if (m_itr_col) {
			cout << setw(4) << N << "|";
		}
		auto it = m_precisions.begin();
		for (const auto& value : values)
		{
			if (value != 0.0 && abs(value) < 0.0001)
				cout << scientific << setprecision(15); else cout << fixed << setprecision(*it);
			cout << setw(*it * 2) << value << "|";
			it++;
		}
		cout << "\n";
	}

	//распечатка строки таблицы
	void printTableRow(const vector<double>& values)
	{
		printTableRow(0, values);
	}
};

//Функция и её производная
class Function
{
	std::function<double(double)> m_func;
	std::function<double(double, int)> m_derivative;
public:
	Function(std::function<double(double)> func, std::function<double(double, int)> derivative = [](double x, int y) { return 0; })
		: m_func(func), m_derivative(derivative)
	{}

	double func(double x) const {
		return m_func(x);
	}

	double derivative(double x, int n = 1) const {
		return m_derivative(x, n);
	}
};

//ДУ 2-ого порядка
//Y''(x) + A(x)Y'(x) - B(x)Y(x) + C(x)×sin(Y(x)) = F(x), B(x) >= 0
class DiffEq
{
public:
	Function m_Ax;
	Function m_Bx;
	Function m_Cx;
	Function m_Fx;

	DiffEq(Function Ax, Function Bx, Function Cx, Function Fx)
		: m_Ax(Ax), m_Bx(Bx), m_Cx(Cx), m_Fx(Fx)
	{}
};

//Решение ДУ 2-ого порядка
class DiffEqSolver
{
public:
	//Таблица результата решения ДУ
	struct Table
	{
		vector<double> m_x;
		vector<double> m_y;
		vector<double> m_z;
		double m_maxErr = 0.0;

		//добавляем новый элемент
		void add(double x, double y, double z) {
			m_x.push_back(x);
			m_y.push_back(y);
			m_z.push_back(z);
		}

		//получить значение искомой функции ДУ в последней точке x = 1.0
		double getY1() {
			return *prev(m_y.end());
		}

		//очистка таблицы
		void clear() {
			m_x.clear();
			m_y.clear();
			m_z.clear();
			m_maxErr = 0.0;
		}
	};

private:
	DiffEq* m_diffEq;
	Table* m_solveTable;
	double m_Y0;
	double m_Z0;
	double m_eps = 0.00001;
	Function* m_Ypr;

public:

	DiffEqSolver(DiffEq* diffEq, double Y0, double Z0, Table* solveTable, Function* Ypr)
		: m_diffEq(diffEq), m_Y0(Y0), m_Z0(Z0), m_solveTable(solveTable), m_Ypr(Ypr)
	{}

	void solve() {
		double x = 0.0;
		double h = 0.1;
		double Y = m_Y0;
		double Z = m_Z0;
		m_solveTable->add(x, Y, Z);
		
		//поиск решения всегда на отрезке [0, 1]
		while (x + h <= 1.0) {
			double nextY1;
			double nextZ1;
			calcNext(x, h, Y, Z, nextY1, nextZ1);

			//используем автоподбор шага h
			double nextY3;
			double nextZ3;
			while(true) {
				//попробуем уточнить, сделав вычисления в промежуточной точке
				double nextY2;
				double nextZ2;
				calcNext(x, h / 2.0, Y, Z, nextY2, nextZ2);
				calcNext(x + h / 2.0, h / 2.0, nextY2, nextZ2, nextY3, nextZ3);
				//если достигнута минимально нужная погрешность
				if (pow(nextY1 - nextY3, 2) + pow(nextZ1 - nextZ3, 2) < m_eps * m_eps) {
					break;
				}
				h /= 2.0;
				nextY1 = nextY2;
				nextZ1 = nextZ2;
			}

			Y = nextY1;
			Z = nextZ1;
			x += h;
			m_solveTable->add(x, Y, Z);

			//вычисляем максимальную погрешность на всей области решения
			double error = abs(m_Ypr->func(x) - Y);
			m_solveTable->m_maxErr = max(m_solveTable->m_maxErr, error);
		}
	}

private:
	//вычисление значения функции в следующей точки
	void calcNext(double x, double h, double Y, double Z, double& nextY, double& nextZ) {
		auto f = [&](double x, double Y, double Z) {
			return -m_diffEq->m_Ax.func(x) * Z +
				m_diffEq->m_Bx.func(x) * Y -
				m_diffEq->m_Cx.func(x) * sin(Y) +
				m_diffEq->m_Fx.func(x);
		};
		double K11 = h * Z;
		double K21 = h * f(x, Y, Z);
		double K12 = h * (Z + K21 / 2.0);
		double K22 = h * f(x + h / 2.0, Y + K11 / 2.0, Z + K21 / 2.0);
		double K13 = h * (Z + K22 / 2.0);
		double K23 = h * f(x + h / 2.0, Y + K12 / 2.0, Z + K22 / 2.0);
		double K14 = h * (Z + K23);
		double K24 = h * f(x + h, Y + K13, Z + K23);

		nextY = Y + (K11 + 2 * K12 + 2 * K13 + K14) / 6.0;
		nextZ = Z + (K21 + 2 * K22 + 2 * K23 + K24) / 6.0;
	}
};

//Метод стрельб
class ShootingMethod
{
	DiffEq* m_diffEq;
	DiffEqSolver::Table* m_solveTable;
	Function* m_Ypr;
	ConsoleTable* m_consoleTable;
	double m_Y0;
	double m_Y1;
	double m_alpha1;
	double m_alpha2;
	double m_eps = 0.0001;
public:
	ShootingMethod(DiffEq* diffEq, double Y0, double Y1, double alpha1, double alpha2, DiffEqSolver::Table* solveTable, Function* Ypr, ConsoleTable* consoleTable)
		: m_diffEq(diffEq), m_Y0(Y0), m_Y1(Y1), m_alpha1(alpha1), m_alpha2(alpha2), m_solveTable(solveTable), m_Ypr(Ypr), m_consoleTable(consoleTable)
	{}

	void solve() {
		int itr = 1;
		double A = m_alpha1;
		double B = m_alpha2;

		//решение для левого конца A
		DiffEqSolver diffEqSolverA(m_diffEq, m_Y0, A, m_solveTable, m_Ypr);
		diffEqSolverA.solve();
		m_consoleTable->printTableRow(itr++, { A, m_solveTable->getY1(), abs(m_Y1 - m_solveTable->getY1()) });
		//решение для правого конца B
		DiffEqSolver diffEqSolverB(m_diffEq, m_Y0, B, m_solveTable, m_Ypr);
		diffEqSolverB.solve();
		m_consoleTable->printTableRow(itr++, { B, m_solveTable->getY1(), abs(m_Y1 - m_solveTable->getY1()) });

		//поиск начального значения для Z(0) методом половинного деления на отрезке [A, B], чтобы удовлетворить Y(1) = Y_1
		while (B - A >= m_eps) {
			auto Z0 = (A + B) / 2.0;
			m_solveTable->clear();
			DiffEqSolver diffEqSolver(m_diffEq, m_Y0, Z0, m_solveTable, m_Ypr);
			diffEqSolver.solve();
			if (m_solveTable->getY1() < m_Y1)
				A = Z0;
			else B = Z0;
			m_consoleTable->printTableRow(itr++, { Z0, m_solveTable->getY1(), m_solveTable->m_maxErr });
		}
	}
};







//2 ЧАСТЬ
int Variant = 0;
# define M_PI 3.14159265358979323846
double chi = 0.2;

double f_(double t, double x) {
	return 0.1 * Variant * sin(M_PI * x) + 0.1 * chi * M_PI * M_PI * t * Variant * sin(M_PI * x);
}
double u(double t, double x) {
	return x + 0.1 * t * Variant * sin(M_PI * x);
}
double phi(double x) {
	return x;
}

//Тау для явной схемы
double tnExpl(double n, double i) { return (1. / n) * (1. / n) / (4. * 0.2) * i; }

//Тау для неявной схемы
double tnImpl(double n, double i) { return i * (1. / n); }

//Явная схема
vector<double> ExplicitScheme(double chi, double a, double b, int n, int maxi)
{
	vector <double> un(n + 1);
	auto h = 1.0 / n;
	auto tau = h * h / (4 * chi);
	for (auto j = 0; j < un.size(); j++)
	{
		un[j] = phi(j * h);
	}

	auto unnew = un;
	for (auto i = 0; i < maxi; i++)
	{
		for (auto j = 1; j < un.size() - 1; j++)
			unnew[j] = un[j] + tau * (chi * (un[j + 1] - 2 * un[j] + un[j - 1]) / (h * h) + f_(tau * i, j * h));
		un = unnew;
	}
	return un;
}

//Текущая погрешнсть 
static double CurDelta(vector<double> layer, double i, double chi, bool isExplicit)
{
	auto n = layer.size() - 1;
	auto h = 1.0 / n;
	auto tau = isExplicit ? h * h / (4 * chi) : h;
	double max = 0;
	for (int j = 0; j <= n; j++)
	{
		if (abs(u(tau * i, j * h) - layer[j]) > max)
			max = abs(u(tau * i, j * h) - layer[j]);
	}
	return max;
}

//Решение методом прогонки
vector <double> Solve(vector<vector<double>> matrix, vector<double> rightPart)
{
	auto n = matrix[0].size();
	vector <double> alpha(n), betta(n);
	alpha[1] = -matrix[0][1] / matrix[0][0];
	betta[1] = rightPart[0] / matrix[0][0];

	for (int i = 1; i < n - 1; i++)
	{
		alpha[i + 1] = -matrix[i][i + 1] / (matrix[i][i] + matrix[i][i - 1] * alpha[i]);
		betta[i + 1] = (rightPart[i] - matrix[i][i - 1] * betta[i]) / (matrix[i][i] + matrix[i][i - 1] * alpha[i]);
	}

	vector<double> x(n);
	x[n - 1] = (rightPart[n - 1] - matrix[n - 1][n - 2] * betta[n - 1]) / (matrix[n - 1][n - 1] + matrix[n - 1][n - 2] * alpha[n - 1]);
	for (int i = n - 2; i >= 0; i--)
		x[i] = alpha[i + 1] * x[i + 1] + betta[i + 1];
	return x;
}

//Неявная схема
vector<double> ImplicitScheme(double chi, double a, double b, int n, int maxi)
{
	auto h = 1.0 / n;
	auto tau = h;
	auto d = (tau * chi) / (h * h);

	vector <double> g(n - 1);
	vector<vector<double>> matrix(n - 1, g);
	for (int i = 0; i < matrix[0].size(); i++)
	{
		if (i - 1 >= 0)
			matrix[i][i - 1] = -d;
		matrix[i][i] = 1 + 2 * d;
		if (i + 1 < matrix[1].size())
			matrix[i][i + 1] = -d;
	}

	vector <double> un(n + 1);
	for (int j = 0; j < un.size(); j++)
	{
		un[j] = phi(j * h);
	}
	un[0] = a;
	un[un.size() - 1] = b;

	vector <double> rightPart(un.size() - 2);
	for (int i = 0; i < maxi; i++)
	{
		for (int j = 0; j < rightPart.size(); j++)
			rightPart[j] = un[j + 1] + tau * f_(tau * (i + 1), h * (j + 1));
		rightPart[0] += d * a;
		rightPart[rightPart.size() - 1] += d * b;
		rightPart = Solve(matrix, rightPart);
		for (int j = 0; j < rightPart.size(); j++)
			un[j + 1] = rightPart[j];
	}
	return un;
}


//Метод конечных разностей
void FiniteDifferenceMethod()
{
	double a(0.), b(1.), Del_T(0);
	vector <double > currentLayer, firstLayer;

	cout << "\nxI = 0.20\nМетод конечных разностей (явная схема)";

	for (int n = 8; n <= 32; n *= 2)
	{
		cout << "\nN = " << n << endl << endl;
		cout << setw(10) << "t" << setw(26) << "delta" << setw(4) << "x:\n";
		Del_T = 0;
		bool last = false;
		for (int i = 1; round(tnExpl(n, i) * 1000.) / 1000. <= 1.001; i++)
		{
			double t = tnExpl(n, i);
			currentLayer = ExplicitScheme(chi, a, b, n, i);
			if (i == 0) firstLayer = currentLayer;
			auto currentDelta = CurDelta(currentLayer, i, chi, true);

			// Output
			cout << setw(7) << setprecision(3) << fixed << t << "|" << setw(22) << setprecision(18);
			if (currentDelta < 0.0001)  cout << setprecision(14) << scientific << currentDelta << setw(4) << "| ";
			else cout << fixed << currentDelta << setw(4) << "| ";
			for (auto i : currentLayer)
				cout << setw(9) << setprecision(5) << fixed << i << "|";
			cout << endl;

			// maxDelta
			if (currentDelta > Del_T)
				Del_T = currentDelta;
			if (round(t * 1000.) >= 1000) break;
		}

		cout << setprecision(15) << fixed << "\nDel_T: " << Del_T << endl << endl;
	}

	cout << "\nМетод конечных разностей (неявная схема)";

	for (int n = 8; n <= 32; n *= 2)
	{
		cout << "\nN = " << n << endl << endl;
		cout << setw(10) << "t" << setw(26) << "delta" << setw(4) << "x:\n";

		Del_T = 0;
		for (int i = 1; tnImpl(n, i) <= 1; i++)
		{
			double t = tnImpl(n, i);
			currentLayer = ImplicitScheme(chi, a, b, n, i);
			if (i == 0) firstLayer = currentLayer;
			auto currentDelta = CurDelta(currentLayer, i, chi, false);

			// Output
			cout << setw(7) << setprecision(3) << fixed << t << "|" << setw(22) << setprecision(18);
			if (currentDelta < 0.0001)  cout << setprecision(14) << scientific << currentDelta << setw(4) << "| ";
			else cout << fixed << currentDelta << setw(4) << "| ";
			for (auto i : currentLayer)
				cout << setw(9) << setprecision(5) << fixed << i << "|";
			cout << endl;

			// MaxDelta
			if (currentDelta > Del_T)
				Del_T = currentDelta;
		}

		cout << setprecision(15) << fixed << "\nDel_T: " << Del_T;
	}
}

int main()
{
	system("chcp 1251");
	
	/*Variant = 25;

	Function Ax([](double x) { return 50.0 * (x + 1.0); });
	Function Bx([](double x) { return x * x + 2.0; });
	Function Cx([](double x) { return x + 1.0; });*/

	Variant = 21;

	Function Ax([](double x) { return 40.0 * (x + 1.0); });
	Function Bx([](double x) { return -x * x + 2.0; });
	Function Cx([](double x) { return x + 1.0; });

	//1 задача
	Function Ypr(
		[&](double x) { return 1 + x + 10 * log(Variant + 1.0) * pow(x, 3) * pow((1 - x), 3); },
		[&](double x, int n) {
			if (n == 1)
				return 1 - 30 * log(Variant + 1.0) * pow(x, 2) * pow((x - 1), 2) * (2 * x - 1);
			if (n == 2)
				return -60 * log(Variant + 1.0) * x * (x - 1) * (5 * pow(x, 2) - 5 * x + 1);
			return 0.0;
		}
	);
	Function Fx([&](double x) { return Ypr.derivative(x, 2) + Ax.func(x) * Ypr.derivative(x, 1) - Bx.func(x) * Ypr.func(x) + Cx.func(x) * sin(Ypr.func(x)); });
	DiffEq diffEq(Ax, Bx, Cx, Fx);

	printf("Метод стрельб\n");
	ConsoleTable consoleTable1({ "z(0)", "y(1)", "Delta" }, { 5, 6, 17 }, true);
	consoleTable1.printTableHeader();

	DiffEqSolver::Table solveTable;
	ShootingMethod shootingMethod(&diffEq, 1.0, 2.0, 0.0, 3.0, &solveTable, &Ypr, &consoleTable1);
	shootingMethod.solve();
	
	ConsoleTable consoleTable2({ "x", "y(x)", "Ypr(x)", "z(x)", "Delta" }, { 5, 5, 6, 5, 17 }, false);
	consoleTable2.printTableHeader();
	for (int i = 0; i < (int)solveTable.m_x.size(); i++) {
		auto x = solveTable.m_x[i];
		auto ypr = Ypr.func(x);
		consoleTable2.printTableRow({ x, solveTable.m_y[i], ypr, solveTable.m_z[i], abs(solveTable.m_y[i] - ypr) });
	}

	//2 задача
	chi = 0.2;
	FiniteDifferenceMethod();
}