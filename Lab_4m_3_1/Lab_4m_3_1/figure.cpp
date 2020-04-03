#include "figure.h"
#include <iostream>

//=======<namespace Matrix_Pro>=======//

MatriX Matrix_Pro::Gauss(MatriX &G_)
{
	MatriX G(G_);
	for (int i = 0; i < G.m; i++)//первый проход
	{
		int l;
		for (l = i; l < G.m; l++)
			if (G.A[i][i] != 0)
				break;
			else
				G.shift(i);
		if (l == G.m)
			continue;

		G.II(i, 1.0 / G.A[i][i]);

		for (int j = i; j < G.m - 1; j++)
			if (G.A[j + 1][i] != 0)
				G.III(j + 1, i, -G.A[j + 1][i]);
	}

	for (int i = 0; i < G.m; i++)//второй проход
	{
		for (int j = 0; j < G.m; j++)
		{
			if (G.A[G.m - i - 1][j] != 0)//
			{
				for (int l = i; l < G.m - 1; l++)
				{
					if (G.A[G.m - l - 2][j] != 0)
						G.III(G.m - l - 2, G.m - i - 1, -G.A[G.m - l - 2][G.m - i - 1]);
				}
				break;
			}
		}
	}

	return G;
}

MatriX_Vector Matrix_Pro::solutionGauss(MatriX_Quad &A, MatriX_Vector &B)
{
	return MatriX_Vector(Gauss(A | B).block_vector(A.size()));//вектор решений
}

MatriX_Quad Matrix_Pro::inverseGauss(MatriX_Quad &A)
{
	if (detGauss(A) == 0)
		std::cout << "ERROR: DETERMINANT = 0\n";

	MatriX_Quad E(A);
	E.initIdentity();

	return Gauss(A | E).block(0, A.size(), A.size(), 2 * A.size()).toMatriX_Quad();
}

double Matrix_Pro::detGauss(MatriX_Quad &A)
{
	double det = 1.0;
	MatriX_Quad G(A);

	for (int i = 0; i < G.m; i++)//первый проход
	{
		int l;
		for (l = i; l < G.m; l++)
			if (G.A[i][i] != 0)
				break;
			else
			{
				det *= pow(-1, G.m - i);//
				G.shift(i);
			}

		if (l == G.m)
			continue;

		det *= G.A[i][i];//
		G.II(i, 1.0 / G.A[i][i]);

		for (int j = i; j < G.m - 1; j++)
			if (G.A[j + 1][i] != 0)
				G.III(j + 1, i, -G.A[j + 1][i]);
	}

	for (int i = 0; i < G.m; i++)
		det *= G.A[i][i];//для очистки совести и не только(нужное)

	return det;
}

MatriX_Vector Matrix_Pro::running(MatriX_Quad &A, MatriX_Vector &B)
{
	int n = A.n;

	MatriX_Vector P(n);
	MatriX_Vector Q(n);
	MatriX_Vector X(n);

	double a = 0;
	double b = A(0, 0);
	double c = A(0, 1);
	double d = B[0];

	P.setElem(0, -c / b);
	Q.setElem(0, d / b);

	for (int i = 1; i < n; i++)
	{
		a = A.A[i][i - 1];
		b = A.A[i][i];
		(i < n - 1) ? (c = A.A[i][i + 1]) : (c = 0);
		d = B.A[i][0];

		P.setElem(i, -c / (b + a*P[i - 1]));
		Q.setElem(i, (d - a*Q[i - 1]) / (b + a*P[i - 1]));
	}

	X.setElem(n - 1, Q[n - 1]);
	for (int i = n - 2; i >= 0; i--)
		X.setElem(i, Q[i] + P[i] * X[i + 1]);

	return X;
}

double Matrix_Pro::formulTrapeze(double(*f)(double), double a, double b, double h)
{
	double S = 0;
	double n = (b - a) / h;

	S = h / 2.0 * (f(a) + f(b));
//	std::cout << a << "\t\t" << f(a) << "\n";
	for (int i = 1; i < n; i++)
	{
		S += h*f(a + h*i);
//		std::cout << a + h*i << "\t\t" << f(a + h*i) << "\n";
	}
//	std::cout << b << "\t\t" << f(b) << "\n";

	return S;
}

double Matrix_Pro::formulSimpson(double(*f)(double), double a, double b, double h)
{
	double S = 0;
	double n = (b - a) / h;

	S = h / 3.0*(f(a) + f(b));
	for (int i = 1; i < n; i++)
	{
		if (i % 2 != 0)
			S += 4.0*h / 3.0*f(a + h*i);
		else if (i < n - 1)
			S += 2.0*h / 3.0*f(a + h*i);
	}

	return S;
}

double Matrix_Pro::formulRunge_Romberg(double z1, double z2, double h1, double h2, int p)
{
	if (h1 != h2)
		return z1 + (z1 - z2) / (pow(h2 / h1, p) - 1.0);
	else
		return z1;
}

double Matrix_Pro::formulRunge(MatriX_Vector &Z, MatriX_Vector &H, int p)
{
	MatriX_Quad D1(Z.size());
	MatriX_Quad D2(Z.size());

	for (int i = 0; i < D1.size(); i++)
	{
		D1.A[i][0] = Z[i];
		D2.A[i][0] = 1;
		for (int j = 1; j < D1.size(); j++)
		{
			D1.A[i][j] = pow(H[i], p + j - 1);
			D2.A[i][j] = pow(H[i], p + j - 1);
		}
	}

	return detGauss(D1) / detGauss(D2);
}

void Matrix_Pro::methodEuler_ODU_I(MatriX_Vector &X, MatriX_Vector &Y, double (*f)(double, double), double a, double b, double y0, double h)
{
	double n = (b - a) / h;
	X = MatriX_Vector(int(n + 1.1));
	Y = MatriX_Vector(int(n + 1.1));

	X.A[0][0] = a;
	Y.A[0][0] = y0;
	for (int i = 0; i < n; i++)
	{
		X.A[i + 1][0] = X[i] + h;
		Y.A[i + 1][0] = Y[i] + h*f(X[i], Y[i]);
	}
}

void Matrix_Pro::methodRunge_Kutta_ODU_I(MatriX_Vector &X, MatriX_Vector &Y, double(*f)(double, double), double a, double b, double y0, double h)
{
	double n = (b - a) / h;
	X = MatriX_Vector(int(n + 1.1));
	Y = MatriX_Vector(int(n + 1.1));

	X.A[0][0] = a;
	Y.A[0][0] = y0;
	for (int i = 0; i < n; i++)
	{
		double k1 = f(X[i], Y[i]);
		double k2 = f(X[i] + h / 2.0, Y[i] + h / 2.0*k1);
		double k3 = f(X[i] + h / 2.0, Y[i] + h / 2.0*k2);
		double k4 = f(X[i] + h, Y[i] + h*k3);

//		std::cout << k1 << "\t" << k2 << "\t" << k3 << "\t" << k4 << "\n";

		X.A[i + 1][0] = X[i] + h;
		Y.A[i + 1][0] = Y[i] + h / 6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
	}
}

void Matrix_Pro::methodEuler_ODU_II(MatriX_Vector &X, MatriX_Vector &Y, MatriX_Vector &Z, 
	double(*f)(double, double, double), double a, double b, double y0, double z0, double h)
{
	double n = (b - a) / h;
	X = MatriX_Vector(int(n + 1.1));
	Y = MatriX_Vector(int(n + 1.1));
	Z = MatriX_Vector(int(n + 1.1));

	X.A[0][0] = a;
	Y.A[0][0] = y0;
	Z.A[0][0] = z0;
	for (int i = 0; i < n; i++)
	{
		X.A[i + 1][0] = X[i] + h;
		Y.A[i + 1][0] = Y[i] + h*Z[i];
		Z.A[i + 1][0] = Z[i] + h*f(X[i], Y[i], Z[i]);
	}
}

void Matrix_Pro::methodRunge_Kutta_ODU_II(MatriX_Vector &X, MatriX_Vector &Y, MatriX_Vector &Z, 
	double(*f)(double, double, double), double a, double b, double y0, double z0, double h)
{
	double n = (b - a) / h;
	X = MatriX_Vector(int(n + 1.1));
	Y = MatriX_Vector(int(n + 1.1));
	Z = MatriX_Vector(int(n + 1.1));

	X.A[0][0] = a;
	Y.A[0][0] = y0;
	Z.A[0][0] = z0;
	for (int i = 0; i < n; i++)
	{
		double ky1 = Z[i];					double kz1 = f(X[i], Y[i], Z[i]);
		double ky2 = Z[i] + h / 2.0*kz1;	double kz2 = f(X[i] + h / 2.0, Y[i] + h / 2.0*ky1, Z[i] + h / 2.0*kz1);
		double ky3 = Z[i] + h / 2.0*kz2;	double kz3 = f(X[i] + h / 2.0, Y[i] + h / 2.0*ky2, Z[i] + h / 2.0*kz2);
		double ky4 = Z[i] + h*kz3;			double kz4 = f(X[i] + h, Y[i] + h*ky3, Z[i] + h*kz3);

		std::cout.width(13);
		std::cout << "\n# " << i << " ky: " << ky1 << " " << ky2 << " " << " " << ky3 << " " << ky3 << "\n";
		std::cout << "# " << i << " kz: " << kz1 << " " << kz2 << " " << " " << kz3 << " " << kz3 << "\n\n";

		X.A[i + 1][0] = X[i] + h;
		Y.A[i + 1][0] = Y[i] + h / 6.0*(ky1 + 2 * ky2 + 2 * ky3 + ky4);
		Z.A[i + 1][0] = Z[i] + h / 6.0*(kz1 + 2 * kz2 + 2 * kz3 + kz4);
	}
}

void Matrix_Pro::methodFinDiff_ODU_II(MatriX_Vector &X, MatriX_Vector &Y, 
	double(*k)(double), double(*l)(double), double(*m)(double), double(*f)(double),
	double r, double s, double t,
	double v, double w, double z,
	double a, double b, double h)
{
	double n = (b - a) / h + 1;
	X = MatriX_Vector(int(n + 0.1));
	
	MatriX_Quad M(X.size());
	MatriX_Vector D(X.size());

	for (int i = 0; i < X.size(); i++)
		X.A[i][0] = a + h*i;

	M.A[0][0] = -r / h + s;
	M.A[0][1] = r / h;
	D.A[0][0] = t;
	for (int i = 1; i < M.size() - 1; i++)
	{
		M.A[i][i - 1] = k(X[i]) / (h*h) - l(X[i]) / (2.0*h);
		M.A[i][i] = -2.0*k(X[i]) / (h*h) + m(X[i]);
		M.A[i][i + 1] = k(X[i]) / (h*h) + l(X[i]) / (2.0*h);
		D.A[i][0] = f(X[i]);
	}
	M.A[M.size() - 1][M.size() - 2] = v / h;
	M.A[M.size() - 1][M.size() - 1] = -v / h - w;
	D.A[D.size() - 1][0] = -z;

//	(M | D).print(5);

	Y = running(M, D);
}

//=======</namespace Matrix_Pro>=======//












