#pragma once
#include <iostream>

namespace TME
{
	//����� ����������� ������� (���������)
	double dichotomy(double (*f)(double), double a_, double b_, double e);

	//����� ���� (�������)
	double chord(double (*f)(double), double a_, double b_, double e);

	//����� ����������� (�������)
	double Newtone(double(*f)(double), double(*f_diff)(double), double(*f_diff2)(double), double a_, double b_, double e);

	//����� ��������
	double iter(double (*f)(double), double (*g)(double), double (*g_diff)(double), double a_, double b_, double e);
}