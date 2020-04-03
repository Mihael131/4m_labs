#include <iostream>
#include <iomanip>
#include <locale>
#include "figure.h"

using namespace std;
using namespace Matrix_Pro;

MatriX_Quad Jacobian1(MatriX_Vector &MX)//для каждой системы индивидуально
{
	MatriX_Quad MQ(2);

	MQ.A[0][0] = 12 * MX[0] * MX[0];	MQ.A[0][1] = -3 * MX[1] * MX[1];
	MQ.A[1][0] = 3 * MX[1];				MQ.A[1][1] = 3 * MX[0] - 2;

	return MQ;
}

MatriX_Vector F1(MatriX_Vector &MX)//для каждой системы индивидуально
{
	MatriX_Vector MF(2);

	MF.A[0][0] = 4 * MX[0] * MX[0] * MX[0] - MX[1] * MX[1] * MX[1] + 19;
	MF.A[1][0] = 3 * MX[0] * MX[1] - 2 * MX[1] - 20;

	return MF;
}

MatriX_Quad Jacobian2(MatriX_Vector &MX)//для каждой системы индивидуально
{
	MatriX_Quad MQ(3);

	MQ.A[0][0] = -2;		MQ.A[0][1] = 5;			MQ.A[0][2] = 4;
	MQ.A[1][0] = 6 * MX[0];	MQ.A[1][1] = 6 * MX[1];	MQ.A[1][2] = 2;
	MQ.A[2][0] = 4 * MX[1];	MQ.A[2][1] = 3 * MX[0];	MQ.A[2][2] = -5;

	return MQ;
}

MatriX_Vector F2(MatriX_Vector &MX)//для каждой системы индивидуально
{
	MatriX_Vector MF(3);

	MF.A[0][0] = -2 * MX[0] + 5 * MX[1] + 4 * MX[2];
	MF.A[1][0] = 3 * MX[0] * MX[0] + 3 * MX[1] * MX[1] + 2 * MX[2] - 5;
	MF.A[2][0] = 4 * MX[0] * MX[1] - 5 * MX[2];

	return MF;
}

int main()
{
	setlocale(LC_ALL, "");
	MatriX_Vector MB;

	cout << setprecision(10);
	cout << "Метод линеризации Ньютона\n\n";
/*
	double A1[2][1] = {
		{2},
		{4}
	};

	double B1[2][1] = {
		{ -2 },
		{ -2 }
	};

	MB.init(&A1[0][0], 2);

	cout << "Уравнение\n";
	cout << "4*x*x*x - y*y*y + 19 = 0\n3*x*y - 2*y - 20 = 0\n\n";
	cout << "Начальное приближение [2,4]\n\nРешение (x,y): ";
	Newtone(F1, Jacobian1, MB, 5).T().print();
	cout << "Проверка (подставляем ответ в исходное уравнение)\n";
	F1(Newtone(F1, Jacobian1, MB, 5)).T().print();

	MB.init(&B1[0][0], 2);

	cout << "Начальное приближение [2,4]\n\nРешение (x,y): ";
	Newtone(F1, Jacobian1, MB, 5).T().print();
	cout << "Проверка (подставляем ответ в исходное уравнение)\n";
	F1(Newtone(F1, Jacobian1, MB, 5)).T().print();
*/
	double A2[3][1] = {
		{ 1 },
		{ 1 },
		{ 1 }
	};

	double B2[3][1] = {
		{ -1 },
		{ -1 },
		{  1 }
	};

	MB.init(&A2[0][0], 3);

	cout << "\nУравнение\n";
	cout << "-2*x + 5*y + 4*z = 0\n";
	cout << "3*x*x + 3*y*y + 2*z - 5 = 0\n";
	cout << "4*x*y - 5*z = 0\n\n";
	cout << "Начальное приближение [1,1,1]\n\nРешение (x,y,z): ";
	Newtone(F2, Jacobian2, MB, 5).T().print();
	cout << "Проверка (подставляем ответ в исходное уравнение)\n";
	F2(Newtone(F2, Jacobian2, MB, 5)).T().print();

	MB.init(&B2[0][0], 3);

	cout << "Начальное приближение [-1,-1,1]\n\nРешение (x,y,z): ";
	Newtone(F2, Jacobian2, MB, 5).T().print();
	cout << "Проверка (подставляем ответ в исходное уравнение)\n";
	F2(Newtone(F2, Jacobian2, MB, 5)).T().print();

//	cout << "\n\nИнтерполяция\n\n";

	double X[5][1] = {
		{ -1 },
		{  0 },
		{  1 },
		{  2 },
		{  3 }
	};

	double Y[5][1] = {
		{ 20 },
		{ 18 },
		{ 20 },
		{ 21 },
		{ 24 }
	};

	MatriX_Vector MX(&X[0][0], 5);
	MatriX_Vector MY(&Y[0][0], 5);
	MatriX_Vector P;
	MatriX_Vector A;
	
/*	cout << "x = ";
	MX.T().print();
	cout << "y = ";
	MY.T().print();

	P = interPolynom_Canon(MX, MY);
	A = interPolynom_Newton(MX.block(1, 5), MY.block(1, 5));
	cout << "Канонический интерполяционный многочлен\n";
	cout << "x = ";	MX.T().print();
	cout << "Коэфф. многочлена "; P.T().print();
	cout << "x = 0.5, y = " << polynom_Canon(P, 0.5) << endl << endl;
	cout << "Интерполяционный многочлен Ньютона\n";
	cout << "x = "; MX.block(1, 5).T().print();
	cout << "Коэфф. многочлена "; A.T().print();
	cout << "x = 0.5, y = " << polynom_Newton(A, MX.block(1, 5), 0.5) << endl << endl;
	cout << "Интерполяционный многочлен Лагранжа\n";
	cout << "x = "; MX.block(0, 4).T().print();
	cout << "x = 0.5, y = " << polynom_Lagrange(MX.block(0, 4), MY.block(0, 4), 0.5) << endl << endl;

	*/

	cout << "\n\nМетод наименьших квадратов\n\n";

	double X1[5][1] = {
		{ -3 },
		{ -1 },
		{ 0 },
		{ 1 },
		{ 2 }
	};

	double Y1[5][1] = {
		{ -7 },
		{ -1 },
		{ 0 },
		{ -1 },
		{ -3 }
	};

	double X11[5][1] = {
		{ -2 },
		{ -1 },
		{ 0 },
		{ 1 },
		{ 2 }
	};

	double Y11[5][1] = {
		{ 3 },
		{ 4 },
		{ 2 },
		{ 1 },
		{ 1 }
	};

	MX.init(&X1[0][0], 5);
	MY.init(&Y1[0][0], 5);

	cout << "x = ";
	MX.T().print();
	cout << "y = ";
	MY.T().print();

	cout << "Квдаратичная апроксимация (a,b,c)\n";
	nApr(MX, MY, 3).print();
	double nev = 0;
	for (int i = 0; i < MY.size(); i++)
		nev += pow(polynom_Canon(nApr(MX, MY, 3), MX[i]) - MY[i], 2);
	cout << "Невязки: " << nev << "\n";

	cout << "\n\nЛинейная аппроксимация (a,b)\n";
	nApr(MX, MY, 2).print();
	nev = 0;
	for (int i = 0; i < MY.size(); i++)
		nev += pow(polynom_Canon(nApr(MX, MY, 2), MX[i]) - MY[i], 2);
	cout << "Невязки: " << nev << "\n\n";

/*	cout << "\n\n\nИнтерполяция сплайнами\n";

	double X2[6][1] = {
		{ -6 },
		{ -3 },
		{ 0 },
		{ 2 },
		{ 5 },
		{ 6 }
	};

	double Y2[6][1] = {
		{ -16 },
		{ 4 },
		{ -7 },
		{ -16 },
		{ -6 },
		{ 8 }
	};

	MX.init(&X2[0][0], 6);
	MY.init(&Y2[0][0], 6);
	
	cout << "    a         b        c        d\n";
	interSpline(MX, MY).T().print();
	
	cout << "Значения в середине каждого отрезка\n";
	for (int i = 0; i < interSpline(MX, MY).sizeColumn(); i++)
		cout << i + 1 << " сплайн, x = " << (MX[i] + MX[i + 1]) / 2.0 << ", y = " << polynom_Canon(interSpline(MX, MY).block_vector(i), (MX[i] + MX[i + 1]) / 2.0) << "\n";
	*/
	system("pause");
}