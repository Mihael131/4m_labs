#include "figure.h"

double TME::dichotomy(double (*f)(double), double a_, double b_, double e)
{
	if (f(a_)*f(b_) > 0)
		std::cout << "ERROR: dichotomy\n\n";

	double a = a_;
	double b = b_;
	double c = (a + b) / 2.0;

	while (abs(b-a) > e)
	{
		if (f(a)*f(c) < 0)
			b = c;
		else
			a = c;

		c = (a + b) / 2.0;

		std::cout << "a=" << a << " f(a)=" << f(a) << " c=" << c << " f(c)=" << f(c) << " b=" << b << " f(b)=" << f(b) << "\n";
	}

	return c;
}

double TME::chord(double (*f)(double), double a_, double b_, double e)
{
	if (f(a_)*f(b_) > 0)
		std::cout << "ERROR: chord\n\n";
	
	double a = a_;
	double b = b_;
	double c = (a*f(b) - b*f(a)) / (f(b) - f(a));
	bool orient = true;//true - right, false - left

	if (f(a)*f(c) < 0)
		orient = false;

	while (abs(f(c)) > e)//можно добавить второе условие (c(n) - c(n-1) < e)
	{
		if (orient)
			a = c;
		else
			b = c;

		c = (a*f(b) - b*f(a)) / (f(b) - f(a));
		std::cout << "a= " << a << " f(a)= " << f(a) << " c= " << c << " f(c)= " << f(c) << " b= " << b << " f(b)= " << f(b) << "\n";
	}

	return c;
}

double TME::Newtone(double(*f)(double), double(*f_diff)(double), double(*f_diff2)(double), double a_, double b_, double e)
{
	if (f(a_)*f(b_) > 0)
		std::cout << "ERROR: Newtone\n\n";

	double a = a_;
	double run;
	double b = b_;

	while (f(a)*f_diff2(a) > 0)//левая граница
	{
		run = a - f(a) / f_diff(a);
		if (abs(a - run) < e)
			return a;
		a = run;
		std::cout << "a= " << a << " f(a)= " << f(a) << "\n";
	}

	while (f(b)*f_diff2(b) > 0)//правая граница
	{
		run = b - f(b) / f_diff(b);
		if (abs(b - run) < e)
			return b;
		b = run;
		std::cout << "b= " << b << " f(b)= " << f(b) << "\n";
	}

	std::cout << "ERROR: Newtone\n\n";
	return 0;
}

double TME::iter(double(*f)(double), double(*g)(double), double(*g_diff)(double), double a_, double b_, double e)
{
	if (f(a_)*f(b_) > 0)
		std::cout << "ERROR: iter\n\n";

	double a = a_;
	double run;
	double b = b_;

	while (g_diff(a) < 1)//левая граница
	{
		run = g(a);
		if (abs(f(a)) < e)
			return a;
		a = run;
		std::cout << "a= " << a << " f(a)= " << f(a) << "\n";
	}

	while (g_diff(b) < 1)//правая граница
	{
		run = g(b);
		if (abs(f(b)) < e)
			return b;
		b = run;
		std::cout << "b= " << b << " f(b)= " << f(b) << "\n";
	}

	std::cout << "ERROR: iter\n\n";

	return 0;
}







