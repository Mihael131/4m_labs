#pragma once
#include "Matrix.h"

namespace Matrix_Pro
{
	MatriX Gauss(MatriX &G);//���������� � ������������ ����
	double detGauss(MatriX_Quad &A);//���������� ������������ ������� ������
	MatriX_Vector solutionGauss(MatriX_Quad &A, MatriX_Vector &B);//������� ���� ������� ������
	MatriX_Quad inverseGauss(MatriX_Quad &A);//���������� �������� ������� ������� ������
	MatriX_Vector running(MatriX_Quad &A, MatriX_Vector &B);
	
	//����� ����������� ������� 
	MatriX_Vector Newtone(MatriX_Vector (*F)(MatriX_Vector&), MatriX_Quad (*J)(MatriX_Vector&), MatriX_Vector &VX, int count_iter);

	//���������������� ��������� (������������)
	MatriX_Vector interPolynom_Canon(MatriX_Vector &X, MatriX_Vector &Y);
	double polynom_Canon(MatriX_Vector &P, double x);

	//���������������� ���������� ��������
	double polynom_Lagrange(MatriX_Vector &X, MatriX_Vector &Y, double x);

	//���������������� ��������� �������
	double f_Newton(MatriX_Vector &X, MatriX_Vector &Y);
	MatriX_Vector interPolynom_Newton(MatriX_Vector &X, MatriX_Vector &Y);
	double polynom_Newton(MatriX_Vector &A, MatriX_Vector &X, double x);

	//����� ���������� ���������
	//�������� �������������
	//������������� n - 1 - ��� �������
	MatriX_Vector nApr(MatriX_Vector &X, MatriX_Vector &Y, int n);

	MatriX interSpline(MatriX_Vector &X, MatriX_Vector &Y);
}

