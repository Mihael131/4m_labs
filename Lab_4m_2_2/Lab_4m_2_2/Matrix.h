#pragma once

const int Max_size = 20;

class MatriX_Quad;
class MatriX_Vector;

class MatriX
{
public:
	bool isQuad();
	bool isVector();
	MatriX_Quad toMatriX_Quad();
	MatriX_Vector toMatriX_Vector();

	//�������������
	MatriX();
	MatriX(MatriX *M_);
	MatriX(int m_, int n_);
	MatriX(double *A_, int m_, int n_);
	void init(double *A_, int m_, int n_);
	void initSize(int m_, int n_);//�������� ������� ������� � ��������
	void initIdentity();//����������� ������� ���������
	void initNull();//��������� ���� ���������
	void initOne();//���������� ������� ���������

	//������������ �������������� ��� ��������
	void I(int i, int j);//������ ������������ �������������� (������ ������ i � j �������)
	void II(int i, double k);//������ ������������ �������������� (��������� i-�� ������ �� �����)
	void III(int i, int j, double k);//������ ������������ �������������� (���������� � i-�� ������ j-�� ����������� �� ����� k)
	void shift(int i);//(��������� ����������) ��������� �� I ������. - ����������� ����� ����� �� ����� � i-�� ������
	
	void falibity();//���������� ������ (��������� ��������� ������ ������� �������� (0.00000001) )
	
	//����
	void setElem(int i, int j, double elem) { A[i][j] = elem; }
	void operator ()(int i, int j, double elem) { A[i][j] = elem; }//���� ��������

	//�����
	int sizeRow() { return m; }//���������� �����
	int sizeColumn() { return n; }//���������� ��������
	double operator ()(int i, int j) { return A[i][j]; }//���������� �������
	void print();//����� �� �����

	//���������� �������� �������� ��� ������
	MatriX operator =(MatriX &M);//"�����" ����������
	MatriX_Quad operator =(MatriX_Quad &MQ);
	MatriX_Vector operator =(MatriX_Vector &MV);
	MatriX operator +(MatriX &M);//��������
	MatriX operator -(MatriX &M);//���������
	MatriX operator *(MatriX &M);//������������
	MatriX operator *(double k);//���������� �� �����
	MatriX operator /(double k);//������� �� �����
	bool operator ==(MatriX &M);//�������� �� ���������
	MatriX operator|(MatriX &M);//������������ (������������) ������

	//�������� ������
	MatriX block(int i1, int i2, int j1, int j2);//���� ������� � i1 ���. (���) �� i2 ������ (�� ���), ���������� �� �����. j1 j2
	MatriX_Vector block_vector(int j);//���������� j �� ������� 
	MatriX T();//���������� ����������������� �������
	double maxElem();//������������ ������� � �������
	double absMaxElem();//������������ �� ������ ������� � �������
	double minElem();//����������� �������
	double absMinElem();//����������� �� ������ �������

	~MatriX() {}

	double A[Max_size][Max_size];//�������
	int m;//���-�� �����
	int n;//���-�� ��������
};

class MatriX_Quad : public MatriX
{
public:
	MatriX_Quad() : MatriX() {}
	MatriX_Quad(int n_) : MatriX(n_, n_) {}
	MatriX_Quad(double *A_, int n_) : MatriX(A_, n_, n_) {}
	MatriX_Quad(MatriX_Quad *MQ_) : MatriX(MQ_) {}
	void init(double *A_, int n_);
	void initSize(int m_);
	MatriX toMatriX();
	MatriX_Quad T();
	int size() { return m; }

	//���������   ����������� ��� ��� ������ �����
	MatriX operator *(MatriX &M);
	MatriX_Quad operator *(MatriX_Quad &MQ);
	MatriX_Vector operator *(MatriX_Vector &MV);
	MatriX_Quad operator *(double k);
	MatriX_Quad operator /(double k);
	MatriX_Quad operator +(MatriX_Quad &MQ);
	MatriX_Quad operator -(MatriX_Quad &MQ);

	//���� ����������� ���� ��� �������� � �������� ������, �� �� �������� ������������ �� ������������� ���������
};

class MatriX_Vector : public MatriX
{
public:
	MatriX_Vector() : MatriX() {}
	MatriX_Vector(int n_) : MatriX(n_, 1) {}
	MatriX_Vector(double *A_, int n_) : MatriX(A_, n_, 1) {}
	MatriX_Vector(MatriX_Vector *MV_) : MatriX(MV_) {}
	void init(double *A_, int n_);
	void initSize(int m_);
	MatriX toMatriX();
	MatriX_Vector block(int i1, int i2);
	int size() { return m; }

	void setElem(int i, double elem) { A[i][0] = elem; }

	double operator[](int i) { return A[i][0]; }
	MatriX_Vector operator *(double k);
	MatriX operator *(MatriX &M);
	MatriX_Vector operator /(double k);
	MatriX_Vector operator +(MatriX_Vector &MV);
	MatriX_Vector operator -(MatriX_Vector &MV);
	void operator ()(int i, double elem) { A[i][0] = elem; }
};
