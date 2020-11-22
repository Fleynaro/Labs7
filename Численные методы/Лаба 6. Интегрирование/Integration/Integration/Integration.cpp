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
		cout << "|" << setw(4) << "N |";
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
		cout << "|" << setw(4) << N << "|";
		auto it = m_precisions.begin();
		for (const auto& value : values)
		{
			cout << setw(*it * 2) << setprecision(*it) << value << "|";
			it++;
		}
		cout << "\n";
	}
};

//Функция и её производная
class Function
{
	std::function<double(double)> m_func;
	std::function<double(double, int)> m_derivative;
public:
	int m_callingCount = 0;

	Function(std::function<double(double)> func, std::function<double(double, int)> derivative)
		: m_func(func), m_derivative(derivative)
	{}

	double func(double x) {
		m_callingCount++;
		return m_func(x);
	}

	double derivative(double x, int n = 1) {
		return m_derivative(x, n);
	}
};

//Отрезок
struct Interval
{
	double m_A;
	double m_B;

	Interval() = default;
};

//Метод трапеций
class TrapezMethod
{
	Function* m_func;
	Interval m_interval;
	bool m_modificiated;
	Table* m_outTable;
	double m_eps = 0.00000001;
public:
	TrapezMethod(Function* func, Interval interval, bool modificiated , Table* outTable)
		: m_func(func), m_interval(interval), m_modificiated(modificiated), m_outTable(outTable)
	{}

	double calculate() {
		double Sh0; //сумма на предпредыдущем шаге
		double Sh1; //сумма на предыдущем шаге
		double Sh2; //сумма на текущем шаге

		int N = 1;
		double h = m_interval.m_B - m_interval.m_A;
		double sum = (m_func->func(m_interval.m_A) + m_func->func(m_interval.m_B)) / 2.0;
		double dS = m_modificiated * (m_func->derivative(m_interval.m_A) - m_func->derivative(m_interval.m_B)) / 12.0;
		Sh0 = h * (sum + dS * h);
		m_outTable->printTableRow(N, { h, Sh0 });

		N = 2;
		h /= 2.0;
		sum += m_func->func((m_interval.m_A + m_interval.m_B) / 2.0);
		Sh1 = h * (sum + dS * h);
		double err = (Sh1 - Sh0) / (4 * (1 + m_modificiated * 3.0) - 1);
		m_outTable->printTableRow(N, { h, Sh1, err });
		
		double k;
		do {
			N *= 2;
			h /= 2.0;
			
			//сумма значений в новых вершинах
			for (int i = 1; i < N; i += 2) {
				sum += m_func->func(m_interval.m_A + i * h);
			}

			Sh2 = h * (sum + dS * h);
			err = (Sh2 - Sh1) / 3.0;
			k = log((Sh2 - Sh0) / (Sh1 - Sh0) - 1) / log(0.5);
			m_outTable->printTableRow(N, { h, Sh2, err, k });

			Sh0 = Sh1;
			Sh1 = Sh2;
		} while (abs(err / Sh2) >= m_eps);

		return Sh2;
	}
};

//Метод симпсона
class SimpsonMethod
{
	Function* m_func;
	Interval m_interval;
	Table* m_outTable;
	double m_eps = 0.00000001;
public:
	SimpsonMethod(Function* func, Interval interval, Table* outTable)
		: m_func(func), m_interval(interval), m_outTable(outTable)
	{}

	double calculate() {
		double Sh0 = -1; //сумма на предпредыдущем шаге
		double Sh1 = -1; //сумма на предыдущем шаге
		double Sh2; //сумма на текущем шаге

		int N = 1;
		double h = m_interval.m_B - m_interval.m_A;
		double sum0 = m_func->func(m_interval.m_A) + m_func->func(m_interval.m_B);

		double err = 100000.0;
		double k;
		do {
			N *= 2;
			h /= 2.0;

			auto sum = sum0;
			for (int i = 1; i < N; i += 2) {
				sum += 4.0 * m_func->func(m_interval.m_A + i * h);
			}

			for (int i = 2; i < N; i += 2) {
				sum += 2.0 * m_func->func(m_interval.m_A + i * h);
			}

			Sh2 = h * sum / 3;

			if (Sh0 != -1 && Sh1 != -1) {
				err = (Sh2 - Sh1) / 15.0;
				k = log((Sh2 - Sh0) / (Sh1 - Sh0) - 1) / log(0.5);
				m_outTable->printTableRow(N / 2, { h, Sh2, err, k });
			}
			else {
				m_outTable->printTableRow(N / 2, { h, Sh2 });
			}

			Sh0 = Sh1;
			Sh1 = Sh2;
		} while (abs(err / Sh2) >= m_eps);

		return Sh2;
	}
};


int main()
{
	system("chcp 1251");

	Interval interval;
	interval.m_A = 1.0;
	interval.m_B = 2.0;

	//функция и её n-ая производная
	Function func(
		[](double x) { return pow(2, x) + 2 - 3 * x; },
		[](double x, int n) { return pow(2, x) * pow(log(2), n) + (n == 1 ? -1 : 0); }
	);

	Table table({ "h", "Integral", "Оценка. погр.", "k" }, { 4, 6, 17, 5 });


	printf("Формула трапеций 0\n");
	table.printTableHeader();
	TrapezMethod trapezMethod(&func, interval, false, &table);
	auto result = trapezMethod.calculate();
	printf("Результат %.15f\nKobr = %i\n", result, func.m_callingCount);
	func.m_callingCount = 0;

	printf("Формула трапеций 1\n");
	table.printTableHeader();
	TrapezMethod trapezMethodMod(&func, interval, true, &table);
	result = trapezMethodMod.calculate();
	printf("Результат %.15f\nKobr = %i\n", result, func.m_callingCount);
	func.m_callingCount = 0;

	printf("Формула Симпсона\n");
	table.printTableHeader();
	SimpsonMethod simpsonMethod(&func, interval, &table);
	result = simpsonMethod.calculate();
	printf("Результат %.15f\nKobr = %i\n", result, func.m_callingCount);
	func.m_callingCount = 0;
}