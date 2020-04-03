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

MatriX_Vector Matrix_Pro::Newtone(MatriX_Vector (*F)(MatriX_Vector&), MatriX_Quad (*J)(MatriX_Vector&), MatriX_Vector &VX, int count_iter)
{
	MatriX_Vector X(VX);
	MatriX_Vector X1(VX);

	for (int i = 0; i < count_iter; i++)
	{
		X1 = X - inverseGauss(J(X))*F(X);
		X = X1;
	}

	return X;
}

MatriX_Vector Matrix_Pro::interPolynom_Canon(MatriX_Vector &X, MatriX_Vector &Y)
{
	MatriX_Quad A(X.m);

	for (int i = 0; i < X.m; i++)
		for (int j = 0; j < X.m; j++)
			A.A[i][j] = pow(X[i], X.m - 1 - j);

	return solutionGauss(A, Y);
}

double Matrix_Pro::polynom_Canon(MatriX_Vector &P, double x)
{
	double y = 0;

	for (int i = 0; i < P.m; i++)
		y += P[i] * pow(x, P.m - 1 - i);

	return y;
}

double Matrix_Pro::polynom_Lagrange(MatriX_Vector &X, MatriX_Vector &Y, double x)
{
	double sum = 0;
	for (int i = 0; i < X.size(); i++)
	{
		double p = 1;
		for (int j = 0; j < X.size(); j++)
			if (i != j)
				p *= (x - X[j]) / (X[i] - X[j]);
		sum += Y[i] * p;
	}
	return sum;
}

double Matrix_Pro::f_Newton(MatriX_Vector &X, MatriX_Vector &Y)
{
	if (X.m > 1)
		return (f_Newton(X.block(1, X.m), Y.block(1, X.m)) - f_Newton(X.block(0, X.m - 1), Y.block(0, X.m - 1))) / (X[X.m - 1] - X[0]);
	else
		return Y[0];
}

MatriX_Vector Matrix_Pro::interPolynom_Newton(MatriX_Vector &X, MatriX_Vector &Y)
{
	MatriX_Vector A(X);

	for (int i = 0; i < A.m; i++)
		A.A[i][0] = f_Newton(X.block(0, i + 1), Y.block(0, i + 1));

	return A;
}

double Matrix_Pro::polynom_Newton(MatriX_Vector &A, MatriX_Vector &X, double x)
{
	double y = 0;
	for (int i = 0; i < A.m; i++)
	{
		double p = 1;
		for (int j = 0; j < i; j++)
			p *= (x - X[j]);

		y += A[i] * p;
	}
	return y;
}


MatriX_Vector Matrix_Pro::nApr(MatriX_Vector &X, MatriX_Vector &Y, int n)
{
	MatriX_Quad A(n);
	MatriX_Vector B(n);

	for (int i = 0; i < A.size(); i++)
	{
		for (int l = 0; l < X.size(); l++)
			B.A[i][0] += Y[l] * pow(X[l], n - 1 - i);

		for (int j = 0; j < A.size(); j++)
		{
			for (int l = 0; l < X.size(); l++)
				A.A[i][j] += pow(X[l], 2*(n - 1) - (i + j));
		}
	}

	return solutionGauss(A, B);
}

MatriX Matrix_Pro::interSpline(MatriX_Vector &X, MatriX_Vector &Y)
{
	MatriX_Quad A(X.size() - 2);
	MatriX_Vector B(X.size() - 2);
	MatriX_Vector H(X.size() - 1);
	MatriX_Vector M1(X.size() - 2);
	MatriX_Vector M(X.size());
	MatriX S(4, X.size() - 1);

	for (int i = 0; i < H.size(); i++)
		H.A[i][0] = X[i + 1] - X[i];

	for (int i = 0; i < A.size(); i++)
	{
		if (i > 0)				A.A[i][i - 1] = H[i]/6.0;
								A.A[i][i] = (H[i] + H[i + 1]) / 3.0;
		if (i < A.size() - 1)	A.A[i][i + 1] = H[i + 1] / 6.0;

		B.A[i][0] = (Y[i + 2] - Y[i + 1]) / H[i + 1] - (Y[i + 1] - Y[i]) / H[i];
	}

	M1 = running(A, B);
	
	for (int i = 0; i < M1.size(); i++)
		M.A[i + 1][0] = M1[i];

	for (int i = 0; i < S.n; i++)
	{
		S.A[0][i] = (M[i + 1] - M[i]) / (6.0 * H[i]);
		S.A[1][i] = (X[i + 1] * M[i] - X[i] * M[i + 1]) / (2.0 * H[i]);
		S.A[2][i] = (X[i] * X[i] * M[i + 1] - X[i + 1] * X[i + 1] * M[i] + 2.0 * (Y[i + 1] - Y[i] + H[i] * H[i] * (M[i] - M[i + 1]) / 6.0)) / (2.0 * H[i]);
		S.A[3][i] = (X[i + 1] * X[i + 1] * X[i + 1] * M[i] - X[i] * X[i] * X[i] * M[i + 1]) / (6.0 * H[i]) +
			X[i + 1] * (Y[i] - M[i] * H[i] * H[i] / 6.0) / H[i] - X[i] * (Y[i + 1] - M[i + 1] * H[i] * H[i] / 6.0) / H[i];
	}

	return S;
}

//=======</namespace Matrix_Pro>=======//