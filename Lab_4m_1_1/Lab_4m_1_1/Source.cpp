#include <iostream>
#include <iomanip>
#include <locale>
#include "figure.h"

using namespace std;
using namespace Matrix_Pro;

int main()
{
	setlocale(LC_ALL, "");

	MatriX_Quad MA;
	MatriX_Vector MB;
	

//	cout << "====================<Задание 1>======================\n\n";
	cout << setprecision(10);

	MatriX Q(13,4);
	for (int i = 0; i < 13; i++)
	{
		Q.setElem(i, 0, i);
		Q.setElem(i, 1, pow(0.99, pow(2, i)));
		Q.setElem(i, 2, pow(2.2, pow(2, -i)));
		Q.setElem(i, 3, i != 0 ? pow(i, -i) : 0);
	}

	Q.print(15);


/*	double A[4][4] = {
		{  -7, -9,11,-3 },
		{  -9, 10,-3,11 },
		{ -10,-10,-1,-8 },
		{   9, -7,-7,-6 }
	};

	double B[4][1] = {
		{ -31 },
		{  74 },
		{  -6 },
		{ -40 }
	};

	MA.init(&A[0][0], 4);
	MB.init(&B[0][0], 4);

	cout << "Объединеннная матрица системы (A|B)\n";
	(MA | MB).print(3);

	cout << "Решение СЛАУ методом Гаусса\n";
	solutionGauss(MA, MB).print(3);

	cout << "Обратная матрица матрицы A (метод Гаусса)\n";
	inverseGauss(MA).print(15);

	cout << "Проверка. Умножаем исходную матрицу на обратную.\n";
	(MA * inverseGauss(MA)).print();

	cout << "Определитель матрицы A (метод Гаусса): " << detGauss(MA) << endl << endl;
	cout << "====================<Задание 2>======================\n\n";
	
*/	double A1[4][4] = {
		{ 34, 10,  7,  1 },
		{  1,-32,  2, -2 },
		{ -8, -9,-35,-11 },
		{  3,-10, 10, 42 }
	};

	double B1[4][1] = {
		{ -171 },
		{ -180 },
		{  64  },
		{ 70  }
	};

	double A11[3][3] = {
		{ 20, 4,  -8 },
		{ -3, 15,  5 },
		{ 6, 3, -18 }
	};

	double B11[3][1] = {
		{ 1 },
		{ -2 },
		{ 3 }
	};
	
	MA.init(&A1[0][0], 4);
	MB.init(&B1[0][0], 4);

//	MA.init(&A11[0][0], 3);
//	MB.init(&B11[0][0], 3);
/*
	cout << "Объединеннная матрица системы (A|B)\n";
	(MA|MB).print(3);

	cout << "Решение СЛАУ методом простых итераций (кол-во итераций 10)\n";
	solutionSimpleIter(MA, MB, 10).print(3);


	cout << "Решение СЛАУ методом Зейделя (кол-во итераций 5)\n";
	solutionSeidel(MA, MB, 5).print(3);

	solutionGauss(MA, MB).print(10);

	cout << "====================<Задание 3>======================\n\n";

	double A2[5][5] = {
		{ 20,  6,  0,  0,  0 },
		{ -5,-12, -5,  0,  0 },
		{  0,  6, 18, -3,  0 },
		{  0,  0, -4, 16, -4 },
		{  0,  0,  0, -5,-18 }
	};

	double B2[5][1] = {
		{  64 },
		{  57 },
		{ -57 },
		{ -60 },
		{  79 }
	};

	MA.init(&A2[0][0], 5);
	MB.init(&B2[0][0], 5);

	cout << "Объединеннная матрица системы (A|B)\n";
	(MA | MB).print(3);

	cout << "Решение СЛАУ методом прогонки\n";
	running(MA, MB).print(3);

	cout << "====================<Задание 4>======================\n\n";

	double A3[4][4] = {
		{ -9, 0, 9,  21 },
		{ 27, 9, 1,  0 },
		{ -1, 8, 33, 4 },
		{ -2, 3, 9,  5 }
	};

	double B3[4][1] = {
		{  21 },
		{  12 },
		{ -11 },
		{  12 }
	};

	MA.init(&A3[0][0], 4);
	MB.init(&B3[0][0], 4);

	cout << "Объединеннная матрица системы (A|B)\n";
	(MA | MB).print(3);

	cout << "Решение СЛАУ с помощью LU разложения" << endl;
	solutionLU(MA, MB).print(15);

	cout << "Обратная матрица матрицы A	(LU)\n";
	inverseLU(MA).print(16);

	cout << "Определитель матрицы A (LU): " << detLU(MA) << endl << endl;

	cout << "====================<Задание 5>======================\n\n";

	double A4[3][3] = {
		{-16,  3, 12 },
		{  3,  4, 18 },
		{ 12, 18, -2 }
	};

	double A44[3][3] = {
		{ 17,  1, 1 },
		{ 1,  17, 2 },
		{ 1,   2, 4 }
	};

	MA.init(&A4[0][0], 3);

	cout << "Исходная матрица\n";
	MA.print();
	
/*	for (int i = 1; i < 11; i++)
	{
		cout << "Итерация " << i << "\n";
		cout << "Спектральный радиус:\t"<< stepenRadVec(MA, i).absMaxElem() << "\n\n";
		cout << "Собственный вектор:\t";	
		(stepenRadVec(MA, i)).T().print();
		cout << "Собственный вектор (н):\t";
		(stepenRadVec(MA, i) / stepenRadVec(MA, i).absMaxElem()).T().print();//<100 большая погрешность
		cout << "\n";
	}
	cout << "Итерация " << 100 << "\n";
	cout << "Спектральный радиус:\t" << stepenRadVec(MA, 100).absMaxElem() << "\n\n";
	cout << "Собственный вектор:\t";
	(stepenRadVec(MA, 100)).T().print();
	cout << "Собственный вектор (н):\t";
	(stepenRadVec(MA, 100) / stepenRadVec(MA, 100).absMaxElem()).T().print();//<100 большая погрешность
	cout << "\n";


	cout << "\n\n";

	MatriX_Vector SL;
	MatriX_Quad SV;

	for (int i = 0; i < 10; i++)
	{
		cout << "Итерация " << i << "\n";
		rotateYakobi(SL, SV, MA, i);
		cout << "Собственные значения (Якоби)\n";
		SL.T().print();

		for (int j = 0; j < SL.size(); j++)
		{
			cout << "СВ (транспонированный) для СЗ " << SL[j] << "\n";
			SV.block_vector(j).T().print();
		}
		cout << "\n\n";
	}


	
//	cout << "\n\n\n";
/*	cout << "====================<Задание 6>======================\n\n";

	double A5[3][3] = {
		{ 4, -2,  -2 },
		{ 8, -6,  -2 },
		{-4,  4, -10 }
	};

	MA.init(&A5[0][0], 3);

	cout << "Исходная матрица\n";
	MA.print();

	SL = methodQR(MA, 15);
	SV = methodInvIter(MA, methodQR(MA, 15), 10);

	cout << "Собственные значение методом QR\n";
	SL.T().print();

	cout << "Метод обратных итераций:\n\n";
	for (int i = 0; i < SL.size(); i++)
	{
		cout << "СВ (транспонированный) для СЗ " << SL[i] << "\n";
		SV.block_vector(i).T().print();
	}
	*/
//	cout << "Проверка\n";
//	for (int i = 0; i < MA.size(); i++)
//	{
//		(MA*SV.block_vector(i)).print();
//		(SV.block_vector(i)*SL[i]).print();
//	}

	system("pause");
}