#pragma once
#include "Matrix.h"

namespace Matrix_Pro
{
	MatriX Gauss(MatriX &G);//���������� � ������������ ����
	double detGauss(MatriX_Quad &A);//���������� ������������ ������� ������
	MatriX_Vector solutionGauss(MatriX_Quad &A, MatriX_Vector &B);//������� ���� ������� ������
	MatriX_Quad inverseGauss(MatriX_Quad &A);//���������� �������� ������� ������� ������
	MatriX_Vector running(MatriX_Quad &A, MatriX_Vector &B);
	
	//��������� �������������� //��������  � ����� ��������� ����������
	double formulTrapeze(double(*f)(double), double a, double b, double h);//������� ��������
	double formulSimpson(double(*f)(double), double a, double b, double h);//������� ��������
	
	double formulRunge_Romberg(double z1, double z2, double h1, double h2, int p);//������� �����-��������
	double formulRunge(MatriX_Vector &Z, MatriX_Vector &H, int p);//������� �����

	//������ ���� ��� ��� I �������
	void methodEuler_ODU_I(MatriX_Vector &X, MatriX_Vector &Y, 
		double(*f)(double,double), double a, double b, double y0, double h);//����� ������
	void methodRunge_Kutta_ODU_I(MatriX_Vector &X, MatriX_Vector &Y, 
		double(*f)(double, double), double a, double b, double y0, double h);//����� ����� �����

	//������ ���� ��� ��� II �������
	void methodEuler_ODU_II(MatriX_Vector &X, MatriX_Vector &Y, MatriX_Vector &Z, 
		double(*f)(double,double,double), double a, double b, double y0, double z0, double h);//����� ������
	void methodRunge_Kutta_ODU_II(MatriX_Vector &X, MatriX_Vector &Y, MatriX_Vector &Z, 
		double(*f)(double, double, double), double a, double b, double y0, double z0, double h);//����� ����� �����

	//������� ������ ��� ��� II �������
	void methodFinDiff_ODU_II(MatriX_Vector &X, MatriX_Vector &Y,
		double (*k)(double), double (*l)(double), double (*m)(double), double (*f)(double),
		double r, double s, double t, 
		double v, double w, double z, 
		double a, double b, double h);//����� �������� ���������
}

