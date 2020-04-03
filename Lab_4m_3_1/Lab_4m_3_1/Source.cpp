#include <iostream>
#include <iomanip>
#include <locale>
#include "figure.h"

const double Pi = 3.141592654;
const double e = 2.718281828;

using namespace std;
using namespace Matrix_Pro;

double f1(double x)
{
	return (5.0 * x - 1.0) / (x*x + 1.0);
}

double F1(double x)
{
	return 5.0 / 2.0*log(x*x + 1.0)/log(e) - atan(x);
}

double f2(double x, double y)
{
	return -y - 2 * x*x + 4;
}

double F2(double x)
{
	return -2 * x*x + 4 * x + pow(e, -x) * 163 / pow(e, 8);
}

double f_3(double x, double y, double z)
{
	return -z - 2 * y + 3 * x + 3;
}

double F3(double x)
{
	double c1 = 1.78517226;
	double c2 = 0.05193930883;
	return c1*pow(e, -x / 2.0)*sin(sqrt(7.0) / 2.0*x) + c2*pow(e, -x / 2.0)*cos(sqrt(7.0) / 2.0*x) + 3.0*x / 2.0 + 3.0 / 4.0;
}

double k(double x)
{
	return 2;
}

double l(double x)
{
	return -1;
}

double m(double x)
{
	return -1;
}

double f(double x)
{
	return x*x - x + 1;
}

double F4(double x)
{
	return (50 + 2.0 / 3.0*e) / (11.0*pow(e, -1.0 / 2.0))*pow(e, -x / 2.0) - 1.0 / 3.0*pow(e, x) - x*x + 3 * x - 8;
}


int main()
{
	setlocale(LC_ALL, "");
	
	cout << "=======<Задание 1>======" << endl << endl;
	cout << setprecision(10);

	MatriX_Vector X;
	MatriX_Vector Y;
	MatriX_Vector Z(3);
	MatriX_Vector H(3);
/*	cout << "Интеграл от f(x) = (5*x - 1)/(x*x + 1), на интервале [7,11]" << endl << endl;
	cout << endl << "Формула трапеций" << endl << endl;
	for (int i = 0; i < 3; i++)
	{
		H.A[i][0] = 1.0 / (pow(2.0, i));
		Z.A[i][0] = formulTrapeze(f1, 7, 11, H[i]);
		cout << "Шаг " << H[i] << endl;
		cout << "Результат " << Z[i] << endl << endl;
	}
	cout << "Уточнение по формуле Рунге" << endl;
	cout << "Результат " << formulRunge(Z, H, 2) << endl << endl;
	cout << "Точное значение интеграла (вычисление через первообразную)" << endl;
	cout << "Результат " << F1(11) - F1(7) << endl << endl;

	cout << endl << "Формула Симпсона" << endl << endl;
	cout << "Шаг 1" << endl;
	cout << "Результат " << formulSimpson(f1, 7, 11, 1) << endl << endl;
	cout << "Шаг 0.5" << endl;
	cout << "Результат " << formulSimpson(f1, 7, 11, 0.5) << endl << endl;
	cout << "Уточнение по формуле Рунге-Ромберга" << endl;
	cout << "Результат " << formulRunge_Romberg(formulSimpson(f1, 7, 11, 1), formulSimpson(f1, 7, 11, 0.5), 1, 0.5, 4) << endl << endl;
	cout << "Точное значение интеграла (вычисление через первообразную)" << endl;
	cout << "Результат " << F1(11) - F1(7) << endl << endl;

	
	cout << "\n\n=======<Задание 2>======" << endl << endl;
	cout << "Решение краевой задачи для ОДУ II - ого порядка\n\n";
	cout << "2*y\'\' - y\' - y = x*x - x + 1\n\n";
	cout << "0 <= x <= 1, шаг 0.125\n\n";
	cout << "Краевые условия:\n-y(0) - 2*y\'(0) = 3\n4*y(1) - 3*y\'(1) = -2\n\n";

	cout << "Метод конечных разностей (x,y)\n";
	methodFinDiff_ODU_II(X, Y, k, l, m, f, -2, -1, 3, -3, 4, -2, 0, 1, 0.125);
	(X | Y).print(5);

	for (int i = 0; i < X.size(); i++)
	{
		cout.width(5);
		cout << X[i] << " " << F4(X[i]) << "\n";
	}
*/
	MatriX_Vector Y4(3);
/*	cout << "\n\n=======<Задание 3>======" << endl << endl;
	cout << "Задача Коши для ОДУ I - ого порядка\n\n";
	cout << "y\' = -y - 2*x*x + 4\n\n";
	cout << "-8 <= x <= -4, y(-8) = 3\n\n";

	
	H.setElem(0, 1);
	H.setElem(1, 0.5);
	H.setElem(2, 0.4);
	cout << "Метод Эйлера (x,y)\n\n";
	for (int i = 0; i < 3; i++)
	{
		cout << "Шаг " << H[i] << "\n";
		methodEuler_ODU_I(X, Y, f2, -8, -4, 3, H[i]);
		Y4.setElem(i, Y[Y.size() - 1]);
		(X | Y).print();
	}
	cout << "По формуле Рунге y(-4) = " << formulRunge(Y4, H, 1) << "\n";
	cout << "Точное значение y(-4) = " << F2(-4) << "\n";

	H.initSize(2);
	Y4.initSize(2);
	H.setElem(0, 1);
	H.setElem(1, 0.5);
	cout << "\nМеtод Рунге-Кутты (x,y)\n\n";
	for (int i = 0; i < 2; i++)
	{
		cout << "Шаг " << H[i] << "\n";
		methodRunge_Kutta_ODU_I(X, Y, f2, -8, -4, 3, H[i]);
		Y4.setElem(i, Y[Y.size() - 1]);
		(X | Y).print();
	}
	cout << "По формуле Рунге-Ромберга y(-4) = " << formulRunge_Romberg(Y4[0],Y4[1],H[0],H[1],4) << "\n";
	cout << "Точное значение y(-4) = " << F2(-4) << "\n\n";
*/
	
	cout << "\n\n=======<Задание 4>======" << endl << endl;
	cout << "Задача Коши для ОДУ II - ого порядка\n\n";
	cout << "y\'\' = -y\' - 2*y + 3*x + 3\n\n";
	cout << "-4 <= x <= -2, y(-4) = 6, y\'(-4) = 5\n\n";

	MatriX_Vector Y2(3);
	H.initSize(3);
	H.setElem(0, 0.5);
	H.setElem(1, 0.25);
	H.setElem(2, 0.2);
	cout << "Метод Эйлера (x,y,z)\n\n";
	for (int i = 0; i < 3; i++)
	{
		cout << "Шаг " << H[i] << "\n";
		methodEuler_ODU_II(X, Y, Z, f_3, -4, -2, 6, 5, H[i]);
		Y2.setElem(i, Y[Y.size() - 1]);
		(X | Y | Z).print(13);
	}
	cout << "По формуле Рунге y(-2) = " << formulRunge(Y2, H, 1) << "\n\n";

	
	H.initSize(2);
	Y4.initSize(2);
	H.setElem(0, 0.5);
	H.setElem(1, 0.25);
	cout << "\nМеtод Рунге-Кутты (x,y,z)\n\n";
	for (int i = 0; i < 2; i++)
	{
		cout << "Шаг " << H[i] << "\n";
		methodRunge_Kutta_ODU_II(X, Y, Z, f_3, -4, -2, 6, 5, H[i]);
		Y2.setElem(i, Y[Y.size() - 1]);
		(X | Y | Z).print(13);
	}
	cout << "По формуле Рунге-Ромберга y(-2) = " << formulRunge_Romberg(Y2[0], Y2[1], H[0], H[1], 4) << "\n\n";

	cout << "Аналитическое решение\n\n";
	for (int i = 0; i < X.size(); i++)
	{
		cout.width(5);
		cout << X[i] << " " << F3(X[i]) << "\n";
	}
	cout << endl;

	/*
	MatriX_Quad AAA(2);
	MatriX_Vector BBB(2);

	AAA.setElem(0, 0, -pow(e, 2)*sin(2 * sqrt(7.0))); 
	AAA.setElem(0, 1, pow(e, 2)*cos(2 * sqrt(7.0)));
	AAA.setElem(1, 0, 0.5*pow(e, 2)*sin(2 * sqrt(7.0)) + pow(e, 2)*cos(2 * sqrt(7.0))*sqrt(7.0) / 2.0);
	AAA.setElem(1, 1, -0.5*pow(e, 2)*cos(2 * sqrt(7.0)) + pow(e, 2)*sin(2 * sqrt(7.0))*sqrt(7.0) / 2.0);
	BBB.setElem(0, 45.0 / 4.0);
	BBB.setElem(1, 3.5);

	solutionGauss(AAA, BBB).print();
//	(AAA*solutionGauss(AAA, BBB)).print();
	*/
//1.78517226
//0.05193930883

	
	system("pause");
}