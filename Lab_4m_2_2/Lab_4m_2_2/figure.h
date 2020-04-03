#pragma once
#include "Matrix.h"

namespace Matrix_Pro
{
	MatriX Gauss(MatriX &G);//привидение к ступенчатому виду
	double detGauss(MatriX_Quad &A);//нахождение орпеделителя методом Гаусса
	MatriX_Vector solutionGauss(MatriX_Quad &A, MatriX_Vector &B);//решение слау методом Гаусса
	MatriX_Quad inverseGauss(MatriX_Quad &A);//нахождение обратной матрицы методом Гаусса
	MatriX_Vector running(MatriX_Quad &A, MatriX_Vector &B);
	
	//метод линеризации Ньютона 
	MatriX_Vector Newtone(MatriX_Vector (*F)(MatriX_Vector&), MatriX_Quad (*J)(MatriX_Vector&), MatriX_Vector &VX, int count_iter);

	//интерполяционный многочлен (канонический)
	MatriX_Vector interPolynom_Canon(MatriX_Vector &X, MatriX_Vector &Y);
	double polynom_Canon(MatriX_Vector &P, double x);

	//интерполяционный многочлена Лагранжа
	double polynom_Lagrange(MatriX_Vector &X, MatriX_Vector &Y, double x);

	//интерполяционный многочлен Ньютона
	double f_Newton(MatriX_Vector &X, MatriX_Vector &Y);
	MatriX_Vector interPolynom_Newton(MatriX_Vector &X, MatriX_Vector &Y);
	double polynom_Newton(MatriX_Vector &A, MatriX_Vector &X, double x);

	//метод наименьших квадратов
	//линейная аппроксимация
	//аппроксимация n - 1 - ной степени
	MatriX_Vector nApr(MatriX_Vector &X, MatriX_Vector &Y, int n);

	MatriX interSpline(MatriX_Vector &X, MatriX_Vector &Y);
}

