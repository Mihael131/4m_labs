#pragma once
#include "Matrix.h"

namespace Matrix_Pro
{
	MatriX Gauss(MatriX &G);//привидение к ступенчатому виду
	double detGauss(MatriX_Quad &A);//нахождение орпеделителя методом Гаусса
	MatriX_Vector solutionGauss(MatriX_Quad &A, MatriX_Vector &B);//решение слау методом Гаусса
	MatriX_Quad inverseGauss(MatriX_Quad &A);//нахождение обратной матрицы методом Гаусса
	MatriX_Vector running(MatriX_Quad &A, MatriX_Vector &B);
	
	//численное интегрирование //доделать  с более грамотным разбиением
	double formulTrapeze(double(*f)(double), double a, double b, double h);//формула трапеций
	double formulSimpson(double(*f)(double), double a, double b, double h);//формула Симпсона
	
	double formulRunge_Romberg(double z1, double z2, double h1, double h2, int p);//формула Рунге-Ромберга
	double formulRunge(MatriX_Vector &Z, MatriX_Vector &H, int p);//формула Рунге

	//задача Коши для ОДУ I порядка
	void methodEuler_ODU_I(MatriX_Vector &X, MatriX_Vector &Y, 
		double(*f)(double,double), double a, double b, double y0, double h);//метод Эйлера
	void methodRunge_Kutta_ODU_I(MatriX_Vector &X, MatriX_Vector &Y, 
		double(*f)(double, double), double a, double b, double y0, double h);//метод Рунге Кутты

	//задача Коши для ОДУ II порядка
	void methodEuler_ODU_II(MatriX_Vector &X, MatriX_Vector &Y, MatriX_Vector &Z, 
		double(*f)(double,double,double), double a, double b, double y0, double z0, double h);//метод Эйлера
	void methodRunge_Kutta_ODU_II(MatriX_Vector &X, MatriX_Vector &Y, MatriX_Vector &Z, 
		double(*f)(double, double, double), double a, double b, double y0, double z0, double h);//метод Рунге Кутты

	//краевая задача для ОДУ II пордяка
	void methodFinDiff_ODU_II(MatriX_Vector &X, MatriX_Vector &Y,
		double (*k)(double), double (*l)(double), double (*m)(double), double (*f)(double),
		double r, double s, double t, 
		double v, double w, double z, 
		double a, double b, double h);//метод конечных разностей
}

