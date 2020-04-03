#include "Matrix.h"
#include <iostream>

//=======<class MatriX>=======//

MatriX::MatriX()
{
	m = 0;
	n = 0;
}

MatriX::MatriX(MatriX *M_)
{
	*this = *M_;
}

MatriX::MatriX(int m_, int n_)
{
	m = m_;
	n = n_;

	initNull();
}

MatriX::MatriX(double *A_, int m_, int n_)
{
	init(A_, m_, n_);
}

void MatriX::init(double *A_, int m_, int n_)
{
	m = m_;
	n = n_;

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = *(A_ + i * n + j);
}

void MatriX::initIdentity()
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if (i == j)
				A[i][j] = 1;
			else
				A[i][j] = 0;
}

void MatriX::initNull()
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
				A[i][j] = 0;
}

void MatriX::initOne()
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = 1;
}

bool MatriX::isQuad()
{
	if (m == n)
		return true;
	else
		return false;
}

bool MatriX::isVector()
{
	if (n == 1)
		return true;
	else
		return false;
}

MatriX_Quad MatriX::toMatriX_Quad()
{
	if (!isQuad())
		std::cout << "ERROR: MatriX_Quad MatriX::toMatriX_Quad()   Запрещенный вызов функции!\n\n";
	
	int size = n;//более безопасная версия
	if (m < n)
		size = m;

	MatriX_Quad MQ0(size);
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			MQ0.A[i][j] = A[i][j];

	return MQ0;
}

MatriX_Vector MatriX::toMatriX_Vector()
{
	if (!isVector())
		std::cout << "ERROR: MatriX_Vector MatriX::toMatriX_Vector()   Запрещенный вызов функции!\n\n";
	
	MatriX_Vector MV0(m);
	for (int i = 0; i < m; i++)
		MV0.A[i][0] = A[i][0];

	return MV0;
}

void MatriX::I(int i, int j)
{
	double str[Max_size];
	for (int l = 0; l < n; l++)
		str[l] = A[m - 1][l];

	for (int l = 0; l < n; l++)
		A[j][l] = A[i][l];
	for (int l = 0; l < n; l++)
		A[i][l] = str[l];

//	std::cout << "I\n";
//	print();
//		det *= -1;
}

void MatriX::II(int i, double k)
{
	for (int j = 0; j < n; j++)
		A[i][j] *= k;

	falibity();

//	std::cout << "II\n";
//	print();
	//	det *= 1 / k;
}

void MatriX::III(int i, int j, double k)
{
	for (int l = 0; l < n; l++)
		A[i][l] += A[j][l] * k;

	falibity();

//	std::cout << "III\n";
//	print();
}

void MatriX::shift(int i)
{
	if (i == m - 1)
		return;

	double str[Max_size];
	for (int l = 0; l < n; l++)
		str[l] = A[m - 1][l];

	for (int l = 0; l < n; l++)//
		A[m - 1][l] = A[i][l];//

	for (int j = i; j < m - 2; j++)
		for (int l = 0; l < n; l++)//
			A[j][l] = A[j + 1][l];//

	for (int l = 0; l < n; l++)//
		A[m - 2][l] = str[l];//

//	std::cout << "shift\n";
//	print();

	//	det *= pow(-1, m - i);
}

void MatriX::falibity()
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if (abs(A[i][j]) < 0.000001)
				A[i][j] = 0;
}

void MatriX::print()
{
	falibity();
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			std::cout << A[i][j] << " ";
		std::cout << "\n";
	}
	std::cout << "\n";
}

MatriX MatriX::block(int i1, int i2, int j1, int j2)
{
	MatriX M0(i2 - i1, j2 - j1);

	for (int i = 0; i < i2 - i1; i++)
		for (int j = 0; j < j2 - j1; j++)
			M0.A[i][j] = A[i + i1][j + j1];

	return M0;
}

MatriX_Vector MatriX::block_vector(int j)
{
	MatriX_Vector MV0(m);

	for (int i = 0; i < m; i++)
		MV0.A[i][0] = A[i][j];

	return MV0;
}

MatriX MatriX::T()
{
	MatriX M0(n, m);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			M0.A[i][j] = A[j][i];

	return M0;
}

double MatriX::maxElem()
{
	double max = A[0][0];

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if (A[i][j] > max)
				max = A[i][j];

	return max;
}

double MatriX::absMaxElem()
{
	double aMax = abs(maxElem());
	double aMin = abs(minElem());
	if (aMax < aMin)
		return aMin;
	else
		return aMax;
}

double MatriX::minElem()
{
	double min = A[0][0];

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if (A[i][j] < min)
				min = A[i][j];

	return min;
}

double MatriX::absMinElem()
{
	double aMax = abs(maxElem());
	double aMin = abs(minElem());
	if (aMax > aMin)
		return aMin;
	else
		return aMax;
}

MatriX MatriX::operator=(MatriX &M)
{
	m = M.m;
	n = M.n;

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = M.A[i][j];

	return *this;
}

MatriX_Quad MatriX::operator=(MatriX_Quad &MQ)
{
	m = MQ.m;
	n = MQ.n;

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = MQ.A[i][j];

	return this->toMatriX_Quad();
}

MatriX_Vector MatriX::operator=(MatriX_Vector &MV)
{
	m = MV.m;
	n = MV.n;

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = MV.A[i][j];

	return this->toMatriX_Vector();
}

MatriX MatriX::operator+(MatriX &M)
{
	MatriX M0(m, n);
	if ((M.m == m) && (M.n == n))
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				M0.A[i][j] = A[i][j] + M.A[i][j];
	else
		std::cout << "ERROR: MatriX MatriX::operator+(MatriX &M)   Запрещенный вызов функции!\n\n";
	return M0;
}

MatriX MatriX::operator-(MatriX &M)
{
	MatriX M0(m, n);
	if ((M.m == m) && (M.n == n))
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				M0.A[i][j] = A[i][j] - M.A[i][j];
	else
		std::cout << "ERROR: MatriX MatriX::operator-(MatriX &M)   Запрещенный вызов функции!\n\n";
	return M0;
}

MatriX MatriX::operator*(MatriX &M)
{
	MatriX M0(m, M.n);

	if (n == M.m)
		for (int i = 0; i < m; i++)
			for (int j = 0; j < M.n; j++)
				for (int i1 = 0; i1 < n; i1++)
					M0.A[i][j] += A[i][i1] * M.A[i1][j];
	else
		std::cout << "ERROR: MatriX MatriX::operator*(MatriX &M)   Запрещенный вызов функции!\n\n";

	return M0;
}

MatriX MatriX::operator*(double k)
{
	MatriX M0(m, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			M0.A[i][j] = A[i][j] * k;
	return M0;
}

MatriX MatriX::operator/(double k)
{
	MatriX M0(m, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			M0.A[i][j] = A[i][j] / k;
	return M0;
}

bool MatriX::operator==(MatriX& M)
{
	if ((M.m == m) && (M.n == n))
	{
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				if (A[i][j] != M.A[i][j])
					return false;
	}
	else
		return false;

	return true;
}

MatriX MatriX::operator|(MatriX& M)
{
	MatriX M0(m,n + M.n);
	if (M.m == m)
	{
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				M0.A[i][j] = A[i][j];
		for (int i = 0; i < M.m; i++)
			for (int j = 0; j < M.n; j++)
				M0.A[i][j + n] = M.A[i][j];
	}
	return M0;
}

//=======</class MatriX>=======//

//=======<class MatriX_Quad>=======//

void MatriX_Quad::init(double *A_, int n_)
{
	MatriX::init(A_, n_, n_);
}

MatriX MatriX_Quad::toMatriX()
{
	return *this;
}

MatriX_Quad MatriX_Quad::T()
{
	return MatriX::T().toMatriX_Quad();
}

MatriX MatriX_Quad::operator *(MatriX &M) 
{ 
	return MatriX::operator*(M); 
}

MatriX_Quad MatriX_Quad::operator *(MatriX_Quad &MQ) 
{
	return MatriX::operator*(MQ).toMatriX_Quad(); 
}

MatriX_Vector MatriX_Quad::operator *(MatriX_Vector &MV)
{
	return MatriX::operator*(MV).toMatriX_Vector();
}

MatriX_Quad MatriX_Quad::operator *(double k)
{
	return MatriX::operator*(k).toMatriX_Quad();
}

MatriX_Quad MatriX_Quad::operator /(double k)
{
	return MatriX::operator/(k).toMatriX_Quad();
}

MatriX_Quad MatriX_Quad::operator +(MatriX_Quad &MQ)
{
	return MatriX::operator+(MQ).toMatriX_Quad();
}

MatriX_Quad MatriX_Quad::operator -(MatriX_Quad &MQ)
{
	return MatriX::operator-(MQ).toMatriX_Quad();
}

//=======</class MatriX_Quad>=======//

//=======<class MatriX_Vector>=======//

void MatriX_Vector::init(double *A_, int n_)
{
	MatriX::init(A_, n_, 1);
}

MatriX MatriX_Vector::toMatriX()
{
	return *this;
}

MatriX_Vector MatriX_Vector::block(int i1, int i2)
{
	return MatriX::block(i1, i2, 0, 1).toMatriX_Vector();
}

MatriX_Vector MatriX_Vector::operator *(double k)
{
	return MatriX::operator*(k).toMatriX_Vector();
}

MatriX_Vector MatriX_Vector::operator /(double k)
{
	return MatriX::operator/(k).toMatriX_Vector();
}

MatriX_Vector MatriX_Vector::operator +(MatriX_Vector &MV)
{
	return MatriX::operator+(MV).toMatriX_Vector();
}

MatriX_Vector MatriX_Vector::operator -(MatriX_Vector &MV)
{
	return MatriX::operator-(MV).toMatriX_Vector();
}

//=======</class MatriX_Vector>=======//
