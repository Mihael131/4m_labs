#pragma once
#include "Matrix.h"

namespace Matrix_Pro
{
	MatriX Gauss(MatriX &G);//привидение к ступенчатому виду
	double detGauss(MatriX_Quad &A);//нахождение орпеделителя методом Гаусса
	MatriX_Vector solutionGauss(MatriX_Quad &A, MatriX_Vector &B);//решение слау методом Гаусса
	MatriX_Quad inverseGauss(MatriX_Quad &A);//нахождение обратной матрицы методом Гаусса
	
	MatriX_Vector solutionSimpleIter(MatriX_Quad &A, MatriX_Vector &B, int count_iter);//решение слау методом простых итераций
	MatriX_Vector solutionSeidel(MatriX_Quad &A, MatriX_Vector &B, int count_iter);//решение слау методом Зейделя

	MatriX_Vector running(MatriX_Quad &A, MatriX_Vector &B);//метод прогонки

	void LU(MatriX_Quad &A, MatriX_Quad &L, MatriX_Quad &U);//инициализация L и U (разложение LU)
	MatriX_Vector solutionLU(MatriX_Quad &A, MatriX_Vector &B);//решение слау с помощью LU разложения
	MatriX_Quad inverseLU(MatriX_Quad &A);//нахождение обратной матрицы с помощью LU разложения 
	double detLU(MatriX_Quad &A);//нахождение определителя с помощью LU разложения

	MatriX_Vector stepenRadVec(MatriX_Quad &A, int count_iter);//степенной метод вычисления спектрального радиуса
	void rotateYakobi(MatriX_Vector &SL, MatriX_Quad &SV, MatriX_Quad &A, int count_iter);//метод вращений Якоби для симметричной матрицы

	double scalar(MatriX_Vector &v1, MatriX_Vector &v2);//скалярное произведение векторов
	MatriX_Vector proj(MatriX_Vector &v1, MatriX_Vector &v2);//проекция вектора v2 колинеарно вектору v1
	void QR(MatriX_Quad &A, MatriX_Quad &Q, MatriX_Quad &R);//QR разложение
	MatriX_Vector methodQR(MatriX_Quad &A, int count_iter);//метод QR для нахождения собственных значений
	MatriX_Quad methodInvIter(MatriX_Quad &A, MatriX_Vector &L, int count_iter);//метод обратных итераций
}

