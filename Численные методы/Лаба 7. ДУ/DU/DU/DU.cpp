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
			cout << setw(*it * 2) << setprecision(*it) << value << "|";
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
	struct Table
	{
		vector<double> m_x;
		vector<double> m_y;
		vector<double> m_z;
		double m_maxErr = 0.0;

		void add(double x, double y, double z) {
			m_x.push_back(x);
			m_y.push_back(y);
			m_z.push_back(z);
		}

		double getY1() {
			return *prev(m_y.end());
		}

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
		
		while (x + h <= 1.0) {
			double nextY1;
			double nextZ1;
			calcNext(x, h, Y, Z, nextY1, nextZ1);
			double nextY3;
			double nextZ3;
			while(true) {
				double nextY2;
				double nextZ2;
				calcNext(x, h / 2.0, Y, Z, nextY2, nextZ2);
				calcNext(x + h / 2.0, h / 2.0, nextY2, nextZ2, nextY3, nextZ3);
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

			double error = abs(m_Ypr->func(x) - Y);
			m_solveTable->m_maxErr = max(m_solveTable->m_maxErr, error);
		}
	}

private:

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

		//A
		DiffEqSolver diffEqSolverA(m_diffEq, m_Y0, A, m_solveTable, m_Ypr);
		diffEqSolverA.solve();
		m_consoleTable->printTableRow(itr++, { A, m_solveTable->getY1(), abs(m_Y1 - m_solveTable->getY1()) });
		//B
		DiffEqSolver diffEqSolverB(m_diffEq, m_Y0, B, m_solveTable, m_Ypr);
		diffEqSolverB.solve();
		m_consoleTable->printTableRow(itr++, { B, m_solveTable->getY1(), abs(m_Y1 - m_solveTable->getY1()) });

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

int main()
{
	system("chcp 1251");
	
	int Variant = 25;
	Function Ax([](double x) { return 50.0 * (x + 1.0); });
	Function Bx([](double x) { return x * x + 2.0; });
	Function Cx([](double x) { return x + 1.0; });
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
}