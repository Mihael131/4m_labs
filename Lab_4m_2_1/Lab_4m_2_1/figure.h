#pragma once
#include <iostream>

namespace TME
{
	//метод половинного деления (дихотомия)
	double dichotomy(double (*f)(double), double a_, double b_, double e);

	//метод хорд (секущих)
	double chord(double (*f)(double), double a_, double b_, double e);

	//метод касательных (Ньютона)
	double Newtone(double(*f)(double), double(*f_diff)(double), double(*f_diff2)(double), double a_, double b_, double e);

	//метод итераций
	double iter(double (*f)(double), double (*g)(double), double (*g_diff)(double), double a_, double b_, double e);
}