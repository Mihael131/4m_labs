#include <iostream>
#include <locale>
#include "figure.h"

using namespace std;
using namespace TME;

double f1(double x)
{
	return x*x*x + 2 * x + 2.8;
}

double f2(double x)
{
	return x*x*x + 8 * x - 12;
}

double f3(double x)
{
	return 0.5*x*x*x + 12 * x - 15;
}

double f3_diff(double x)
{
	return 1.5*x*x + 12;
}

double f3_diff2(double x)
{
	return 3*x;
}

double f4(double x)
{
	return 3 * x*x*x + x + 15;
}

double g4(double x)
{
	return -pow((x + 15.0) / 3.0, 1.0 / 3.0);
}

double g4_diff(double x)
{
	return -1.0 / 9.0*pow((x + 15.0) / 3.0, -2.0 / 3.0);
}

int main()
{
	setlocale(LC_ALL, "");

	cout << "Метод половинного деления (дихотомия), x*x*x + 2*x + 2.8 = 0, a = -1, b = 0\n";
	cout << "x = " << dichotomy(f1, -1, 0, 0.0001) << endl << endl;

	cout << "Метод секущих (хорд), x*x*x + 8 * x - 12 = 0, a = 1, b = 2\n";
	cout << "x = " << chord(f2, 1, 2, 0.0001) << endl << endl;

	cout << "Метод касательных (Ньютона), 0.5*x*x*x + 12 * x - 15, a = 1, b = 2\n";
	cout << "x = " << Newtone(f3, f3_diff, f3_diff2, 1, 2, 0.0001) << endl << endl;

	cout << "Метод итераций, 3 * x*x*x + x + 15, a = -2, b = -1\n";
	cout << "x = " << iter(f4, g4, g4_diff, -2, -1, 0.0001) << endl << endl;

	system("pause");
	return 0;
}