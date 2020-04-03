#pragma once

const int Max_size = 20;

class MatriX_Quad;
class MatriX_Vector;

class MatriX
{
public:
	bool isQuad();
	bool isVector();
	MatriX_Quad toMatriX_Quad();//костыль, по возможности переделать
	MatriX_Vector toMatriX_Vector();

	//инициализация
	MatriX();
	MatriX(MatriX *M_);
	MatriX(int m_, int n_);
	MatriX(double *A_, int m_, int n_);
	void init(double *A_, int m_, int n_);
	void initIdentity();//присваивает матрице единичную
	void initNull();//обнуление всех элементов
	void initOne();//заполнение матрицы единицами

	//элементарные преобразования
	void I(int i, int j);//первое элементарное преобразование (меняем строки i и j местами)
	void II(int i, double k);//второе элементарное преобразование (домножаем i-ую строку на число)
	void III(int i, int j, double k);//третье элементарное преобразование (прибавляем к i-ой строке j-ую домноженную на число k)
	void shift(int i);//(небольшое расширение) следствие из I преобр. - циклический сдвиг строк на чиная с i-ой строки
	void falibity();//фильтрация мусора (обнуление элементов меньше заданой величины (0.000001) )
	
	//вывод
	void print();//вывод на экран

	//перегрузка основных операций для матриц
	MatriX operator =(MatriX &M);//"умное" присваение
	MatriX_Quad operator =(MatriX_Quad &MQ);
	MatriX_Vector operator =(MatriX_Vector &MV);
	MatriX operator +(MatriX &M);//сложение
	MatriX operator -(MatriX &M);//вычитание
	MatriX operator *(MatriX &M);//перемножение
	MatriX operator *(double k);//домножение на число
	MatriX operator /(double k);//деление на число
	bool operator ==(MatriX &M);//проверка на равенство
	MatriX operator|(MatriX &M);//конкатенация (объединнение) матриц

	//полезные плюшки
	MatriX block(int i1, int i2, int j1, int j2);//блок матрицы с i1 стр. (вкл) по i2 строку (не вкл), аналогично со столб. j1 j2
	MatriX_Vector block_vector(int j);//возвращает j ый столбец 
	MatriX T();//возвращает транспонированную матрицу
	double maxElem();//максимальный элемент в матрице
	double absMaxElem();//максимальный по модулю элемент в матрице
	double minElem();//минимальный элемент
	double absMinElem();//минимальный по модулю элемент

	~MatriX() {}

	double A[Max_size][Max_size];//матрица
	int m;//кол-во строк
	int n;//кол-во столбцов
};

class MatriX_Quad : public MatriX
{
public:
	MatriX_Quad() : MatriX() {}
	MatriX_Quad(int n_) : MatriX(n_, n_) {}
	MatriX_Quad(double *A_, int n_) : MatriX(A_, n_, n_) {}
	MatriX_Quad(MatriX_Quad *MQ_) : MatriX(MQ_) {}
	void init(double *A_, int n_);
	MatriX toMatriX();
	MatriX_Quad T();
	int size() { return m; }
	
	//операторы   перегрузить все что только можно
	MatriX operator *(MatriX &M);
	MatriX_Quad operator *(MatriX_Quad &MQ);
	MatriX_Vector operator *(MatriX_Vector &MV);
	MatriX_Quad operator *(double k);
	MatriX_Quad operator /(double k);
	MatriX_Quad operator +(MatriX_Quad &MQ);
	MatriX_Quad operator -(MatriX_Quad &MQ);

	//если перегрузить хоть раз оператор в дочернем классе, то он начинает отказываться от родительского оператора сука
};

class MatriX_Vector : public MatriX
{
public:
	MatriX_Vector() : MatriX() {}
	MatriX_Vector(int n_) : MatriX(n_, 1) {}
	MatriX_Vector(double *A_, int n_) : MatriX(A_, n_, 1) {}
	MatriX_Vector(MatriX_Vector *MV_) : MatriX(MV_) {}
	void init(double *A_, int n_);
	MatriX toMatriX();
	MatriX_Vector block(int i1, int i2);
	int size() { return m; }

	double operator[](int i) { return A[i][0]; }
	MatriX_Vector operator *(double k);
	MatriX_Vector operator /(double k);
	MatriX_Vector operator +(MatriX_Vector &MV);
	MatriX_Vector operator -(MatriX_Vector &MV);
};
