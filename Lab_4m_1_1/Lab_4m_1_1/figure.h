#pragma once
#include "Matrix.h"

namespace Matrix_Pro
{
	MatriX Gauss(MatriX &G);//���������� � ������������ ����
	double detGauss(MatriX_Quad &A);//���������� ������������ ������� ������
	MatriX_Vector solutionGauss(MatriX_Quad &A, MatriX_Vector &B);//������� ���� ������� ������
	MatriX_Quad inverseGauss(MatriX_Quad &A);//���������� �������� ������� ������� ������
	
	MatriX_Vector solutionSimpleIter(MatriX_Quad &A, MatriX_Vector &B, int count_iter);//������� ���� ������� ������� ��������
	MatriX_Vector solutionSeidel(MatriX_Quad &A, MatriX_Vector &B, int count_iter);//������� ���� ������� �������

	MatriX_Vector running(MatriX_Quad &A, MatriX_Vector &B);//����� ��������

	void LU(MatriX_Quad &A, MatriX_Quad &L, MatriX_Quad &U);//������������� L � U (���������� LU)
	MatriX_Vector solutionLU(MatriX_Quad &A, MatriX_Vector &B);//������� ���� � ������� LU ����������
	MatriX_Quad inverseLU(MatriX_Quad &A);//���������� �������� ������� � ������� LU ���������� 
	double detLU(MatriX_Quad &A);//���������� ������������ � ������� LU ����������

	MatriX_Vector stepenRadVec(MatriX_Quad &A, int count_iter);//��������� ����� ���������� ������������� �������
	void rotateYakobi(MatriX_Vector &SL, MatriX_Quad &SV, MatriX_Quad &A, int count_iter);//����� �������� ����� ��� ������������ �������

	double scalar(MatriX_Vector &v1, MatriX_Vector &v2);//��������� ������������ ��������
	MatriX_Vector proj(MatriX_Vector &v1, MatriX_Vector &v2);//�������� ������� v2 ���������� ������� v1
	void QR(MatriX_Quad &A, MatriX_Quad &Q, MatriX_Quad &R);//QR ����������
	MatriX_Vector methodQR(MatriX_Quad &A, int count_iter);//����� QR ��� ���������� ����������� ��������
	MatriX_Quad methodInvIter(MatriX_Quad &A, MatriX_Vector &L, int count_iter);//����� �������� ��������
}

