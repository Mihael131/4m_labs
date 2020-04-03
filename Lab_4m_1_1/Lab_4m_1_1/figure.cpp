#include "figure.h"
#include <iostream>

//=======<namespace Matrix_Pro>=======//

MatriX Matrix_Pro::Gauss(MatriX &G_)
{
	MatriX G(G_);
	for (int i = 0; i < G.m; i++)//первый проход
	{
		int l;
		for (l = i; l < G.m; l++)//ищем ненулевой элемент
			if (G.A[i][i] != 0)//если нашли, выходим из цикла
				break;
			else
				G.shift(i);//если нет, делаем сдвиг
		if (l == G.m)//если прокрутили все строки и не нашли ненулевого, переходим к след столбцу
			continue;

		G.II(i, 1.0 / G.A[i][i]);//делим элемент на себя

		for (int j = i; j < G.m - 1; j++)//обнуляем все что под ним
			if (G.A[j + 1][i] != 0)
				G.III(j + 1, i, -G.A[j + 1][i]);
	}

	for (int i = 0; i < G.m; i++)//второй проход обнуляем элементы выше главной диагонали
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

MatriX_Vector Matrix_Pro::solutionSimpleIter(MatriX_Quad &A, MatriX_Vector &B, int count_iter)
{
	MatriX_Quad M(A);
	for (int i = 0; i < M.size(); i++)
	{
		M.II(i, -1.0 / A(i, i));
		M.setElem(i, i, 0);
	}

	MatriX_Vector b(B);
	for (int i = 0; i < b.size(); i++)
		b.II(i, 1.0 / A(i, i));

	MatriX_Vector x1(b);
	MatriX_Vector x2(b);
	x1.initNull();//лиbо b

	MatriX_Vector N1(M.size());//нормы
	MatriX_Vector N2(M.size());
	for (int i = 0; i < M.size(); i++)
		for (int j = 0; j < M.size(); j++)
		{
			N1.setElem(i, N1[i] + M(i, j));// Теперь норм
			N2.setElem(i, N2[i] + M(j, i));
		}

	std::cout << "Норма ||B||1: " << abs(N1.absMaxElem()) << " < 1\n";
	std::cout << "Норма ||B||2: " << abs(N2.absMaxElem()) << " < 1\n\n";

	for (int i = 0; i < count_iter; i++)
	{
		x2 = M*x1 + b;
		x1 = x2;

		std::cout << i + 1 << "-ая итерация x = ";
		x1.T().print();
	}

	return x2;
}

MatriX_Vector Matrix_Pro::solutionSeidel(MatriX_Quad &A, MatriX_Vector &B, int count_iter)
{
	MatriX_Quad B1(A);
	MatriX_Quad B2(A);
	for (int i = 0; i < A.size(); i++)
	{
		B1.II(i, -1.0 / A(i,i));
		B2.II(i, -1.0 / A(i,i));

		for (int j = 0; j < A.size(); j++)
		{
			if (i == j)
				B1.setElem(i, j, 0);

			if (i > j)
				B1.setElem(i, j, 0);
			else
				B2.setElem(i, j, 0);
		}
	}

	MatriX_Vector b(B);
	for (int i = 0; i < b.size(); i++)
		b.II(i, 1.0 / A(i, i));

	MatriX_Vector x1(b);
	MatriX_Vector x2(b);
	x1.initNull();

	MatriX_Vector N1(B1.size());//нормы
	MatriX_Vector N2(B1.size());
	for (int i = 0; i < B1.size(); i++)
		for (int j = 0; j < B1.size(); j++)
		{
			N1.setElem(i, N1[i] + (B1 + B2)(i, j));// Теперь норм
			N2.setElem(i, N2[i] + (B1 + B2)(j, i));
		}

	std::cout << "Норма ||B||1: " << abs(N1.absMaxElem()) << " < 1\n";
	std::cout << "Норма ||B||2: " << abs(N2.absMaxElem()) << " < 1\n\n";


	for (int i = 0; i < count_iter; i++)
	{
		MatriX_Vector s(b + B1*x1);
		MatriX_Quad M(B2.size());
		M.initIdentity();
		M = M - B2;

		for (int k = 0; k < x2.size(); k++)//решение прямым ходом сиситемы (E-B2)*x2 = s
		{
			double sum = 0;
			for (int j = 0; j < k; j++)
				sum += M(k, j)*x2[j];
			x2.setElem(k, s[k] - sum);
		}

		std::cout << i + 1 << "-ая итерация:  ";
		x2.T().print(12);
		
		x1 = x2;
	}

	return x2;
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
		X.setElem(i, Q[i] + P[i]*X[i + 1]);

	return X;
}

void Matrix_Pro::LU(MatriX_Quad &A, MatriX_Quad &L, MatriX_Quad &U)
{
	for (int i = 0; i < A.n; i++)
		for (int j = i; j < A.n; j++)
		{
			double sumU = 0;
			double sumL = 0;
			for (int k = 0; k <= i - 1; k++)
			{
				sumU += L(i, k) * U(k, j);
				sumL += L(j, k)* U(k, i);;
			}
			U.setElem(i, j, A(i, j) - sumU);
			L.setElem(j, i, 1.0 / U(i, i)*(A(j, i) - sumL));
		}
}

MatriX_Vector Matrix_Pro::solutionLU(MatriX_Quad &A, MatriX_Vector &B)
{
	MatriX_Quad L(A.size());
	MatriX_Quad U(A.size());
	
	MatriX_Vector X(B);
	MatriX_Vector Y(B);

	LU(A, L, U);

	std::cout << "\nМатрица L\n";
	L.print(13);

	std::cout << "Матрица U\n";
	U.print(13);

	std::cout << "Вектор решений\n";

	//L*U*X = MB
	for (int i = 0; i < L.size(); i++)//L*Y = MB - прямой ход
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
			sum += L(i,j) * Y[j];
		Y.setElem(i, (B[i] - sum) / L(i,i));
	}

	for (int i = U.size() - 1; i >= 0; i--)//Y = U*X  - обратный ход
	{
		double sum = 0;
		for (int j = i + 1; j < U.size(); j++)
			sum += U(i, j)* X[j];
		X.setElem(i, (Y[i] - sum) / U(i, i));
	}

	return X;
}

MatriX_Quad Matrix_Pro::inverseLU(MatriX_Quad &A)
{
	MatriX_Quad L(A.size());
	MatriX_Quad U(A.size());

	LU(A, L, U);

	MatriX_Quad Y(A);
	MatriX_Quad X(A);
	MatriX_Quad E(A);
	E.initIdentity();

	for (int k = 0; k < A.size(); k++)
	{
		for (int i = 0; i < L.size(); i++)//L*Y = MB - прямой ход
		{
			double sum = 0;
			for (int j = 0; j < i; j++)
				sum += L(i, j)*Y(j, k);
			Y.setElem(i, k, (E(i, k) - sum) / L(i, i));
		}

		for (int i = U.size() - 1; i >= 0; i--)//Y = U*X  - обратный ход
		{
			double sum = 0;
			for (int j = i + 1; j < U.size(); j++)
				sum += U(i, j)*X(j, k);
			X.setElem(i, k, (Y(i, k) - sum) / U(i, i));
		}
	}

	return X;
}

double Matrix_Pro::detLU(MatriX_Quad &A)
{
	MatriX_Quad L(A.size());
	MatriX_Quad U(A.size());

	LU(A, L, U);

	double det = 1;
	for (int i = 0; i < A.n; i++)
		det *= L(i, i) * U(i, i);

	return det;
}

MatriX_Vector Matrix_Pro::stepenRadVec(MatriX_Quad &A, int count_iter)
{
	MatriX_Vector V(A.size());
	MatriX_Vector W(A.size());
	W.initOne();

	for (int i = 0; i < count_iter; i++)
	{
		V = A*W;
		W = V / V.absMaxElem();
	}

	return V;//не нормированный
}

void Matrix_Pro::rotateYakobi(MatriX_Vector &SL, MatriX_Quad &SV, MatriX_Quad &A, int count_iter)
{
	MatriX_Quad A1(A);
	MatriX_Quad H(A);
	MatriX_Quad H1(A);
	H.initIdentity();

	for (int l = 0; l < count_iter; l++)
	{
		//1
		int k = 0;
		int m = 1;
		double max = abs(A1(k, m));

		for (int i = 0; i < A1.size(); i++)
			for (int j = i + 1; j < A1.size(); j++)
				if (abs(A1(i, j)) > max)
				{
					k = i;
					m = j;
					max = abs(A1(k, m));
				}
		
		//2
		double f;
		if (abs(A1(k, k) - A1(m, m)) >= 0.0000001)
			f = 0.5*atan(2.0*(A1(k, m) / (A1(k, k) - A1(m, m))));//abs мб
		else if (A1(k, m) < 0)
				f = -atan(1.0);
			 else
				f = atan(1.0);

//		std::cout << "Угол " << f << "\n";

		//3
		H1.initIdentity();
		H1.setElem(k, k, cos(f));	H1.setElem(k, m, -sin(f));
		H1.setElem(m, k, sin(f));	H1.setElem(m, m, cos(f));

		//4
		A1 = H1.T()*A1*H1;
		H = H*H1;
	}
	
	MatriX_Vector L(A1.size());
	
	H = H.T();
	for (int i = 0; i < L.size(); i++)
	{
		L.setElem(i, A1(i, i));
		H.II(i, 1.0 / H.T().block_vector(i).absMaxElem());//	
//		H.II(i, L[i]);//теперь здеся и собст. значения и вектора
	}
	H = H.T();
	
	SL = L;
	SV = H;

//	L.T().print();
//	H.print();

//	for (int i = 0; i < H.size(); i++)
//	{
//		(A*H.block_vector(i)).print();
//		(H.block_vector(i)*L[i]).print();
//	}
}

double Matrix_Pro::scalar(MatriX_Vector &v1, MatriX_Vector &v2)
{
	if (v1.size() != v2.size())
		std::cout << "ERROR: double Matrix_Pro::scalar(MatriX_Vector *v1, MatriX_Vector *v2)\n\n";

	double sum = 0;
	for (int i = 0; i < v1.size(); i++)
		sum += v1[i] * v2[i];

	return sum;
}

MatriX_Vector Matrix_Pro::proj(MatriX_Vector &v1, MatriX_Vector &v2)
{
	return v1*(scalar(v1, v2) / scalar(v1, v1));
}

void Matrix_Pro::QR(MatriX_Quad &A, MatriX_Quad &Q, MatriX_Quad &R)//работает!
{
	Q = A;
	for (int i = 0; i < A.size(); i++)
	{
		MatriX_Vector sum(A.size());
		for (int j = 0; j < i; j++)
			sum = sum + proj(Q.block_vector(j), A.block_vector(i));

		for (int j = 0; j < A.size(); j++)
			Q(j, i, A(j, i) - sum[j]);
	}

	Q = Q.T();
	for (int i = 0; i < Q.size(); i++)
		if (scalar(Q.T().block_vector(i), Q.T().block_vector(i)) != 0)
			Q.II(i, 1.0 / sqrt(abs(scalar(Q.T().block_vector(i), Q.T().block_vector(i)))));
	
	Q = Q.T();
	R = Q.T()*A;
}

MatriX_Vector Matrix_Pro::methodQR(MatriX_Quad &A, int count_iter)
{
	MatriX_Quad A1(A);
	MatriX_Quad Q;
	MatriX_Quad R;

	for (int i = 0; i < count_iter; i++)
	{
		QR(A1, Q, R);
		A1 = R*Q;
	}

	MatriX_Vector L(A1.size());
	for (int i = 0; i < L.size(); i++)
		L.setElem(i, A1(i, i));

	return L;
}

MatriX_Quad Matrix_Pro::methodInvIter(MatriX_Quad &A, MatriX_Vector &L, int count_iter)
{
	MatriX_Vector X;
	MatriX_Quad V(A);
	MatriX_Quad E(A);
	E.initIdentity();

//	for (int i = 0; i < L.size(); i++)
//		std::cout << detGauss(A - E*L[i]) << "\n";

	for (int i = 0; i < L.size(); i++)
	{
		MatriX_Vector X1(L.size());
		X1.initOne();
		for (int j = 0; j < count_iter; j++)
		{
			X = solutionGauss(A - E*L[i], X1);
			X = X / abs(X.absMaxElem());//нормировка!
			X1 = X;
		}	
		for (int j = 0; j < L.size(); j++)
			V.setElem(j, i, X[j]);
	}

	return V;
}

//=======</namespace Matrix_Pro>=======//