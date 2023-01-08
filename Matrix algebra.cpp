#include "Matrix algebra.h"

// ����� ������� ��������������� ����� ������: ������ ����� (rows), ������ �������� (columns) � ��������� �������� ��������� (content)

int Matrix::getrows()
{
	return rows;
}

int Matrix::getcolumns()
{
	return columns;
}

void Matrix::setcontent(int m, int n) // �� ������� content -- ������� ������� ���������������� �������
{
	rows = m;
	columns = n;
	vector<vector<double>> Theta;
	Theta.resize(m);
	for (int i = 0; i < m; i++)
	{
		Theta[i].resize(n);
		fill(Theta[i].begin(), Theta[i].end(), 0.0);
	}
	content = Theta;
}

//�����, ����������� ��������� ��������� ������ ������ ������
void Matrix::setelement(int i, int j, double a)
{
	if (i < 1 || j < 1 || i > rows || j > columns)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		content[i-1][j-1] = a;
	}
}

//�����, ����������� �������� �� ������ ������
double Matrix::getelement(int i, int j)
{
	if (i < 1 || j < 1 || i > rows || j > columns)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		return content[i-1][j-1];
	}
}

//��������� ����������� ������ ������� 3x3
Matrix::Matrix()
{
	rows = 3;
	columns = 3;
	setcontent(3,3);
}
	
//�����������, �������� ������� ������� m x n
Matrix::Matrix(int m, int n)
{
	rows = m;
	columns = n;
	setcontent(m, n);
}

Matrix::Matrix(int m, int n, double k)
{
	rows = m;
	columns = n;
	vector<vector<double>> E;
	E.resize(m);
	for (int i = 0; i < m; i++)
	{
		E[i].resize(n);
		fill(E[i].begin(), E[i].end(), 0.0);
		if (i < n)
		{
			E[i][i] = k;
		}
	}
	content = E;
}

//�����������, �������� ������� m x n � ����������� � ���������� ������� ������� C
Matrix::Matrix(vector<vector<double>>& C)
{
	int m = C.size();
	int n = C[0].size();
	setcontent(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			content[i][j] = C[i][j];
		}
	}
	cout << "Matrix created!" << endl;
}

//�������� ������� � ������
void Vector::setelement(int i, double a)
{
	if (i < 1 || i > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		content[i-1][0] = a;
	}
}

//�������� ������� �������
double Vector::getelement(int i)
{
	if (i < 1 || i > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		return content[i-1][0];
	}
}

double Vector::norm()
{
	double s = 0;
	for (int i = 0; i < rows; i++)
	{
		s += content[i][0] * content[i][0];
	}
	return sqrt(s);
}

//��������� ����������� ������� -- ������� ������� 3 x 1
Vector::Vector()
{
	rows = 3;
	columns = 1;
	setcontent(3, 1);
}

//������� ������ �� m �������
Vector::Vector(int m)
{
	rows = m;
	columns = 1;
	setcontent(m, 1);
}

Vector::Vector(const Vector& copia)
{
	rows = copia.rows;
	content = copia.content;
}

//������ ������� m, ����������� ���������� ������� C
Vector::Vector(vector<double>& C)
{
	int m = C.size();
	columns = 1;
	setcontent(m,1);
	for (int i = 0; i < m; i++)
	{
		content.at(i).at(0) = C.at(i);
	}
}

//�������� ������������ ��� ��������
Vector& Vector::operator= (const Vector& v)
{
	if (rows != v.rows)
	{
		throw "Error.";
	}
	for (int i = 0; i < v.rows; i++)
	{
		content[i].at(0) = v.content[i].at(0);
	}
	return *this;
}

void Matrix::permutation(int i, int j)
{
	if (i < 1 || i > rows || j  < 1 || j > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		vector <double> buf = content[i - 1];
		content[i - 1] = content[j - 1];
		content[j - 1] = buf;
	}
}

void Matrix::scalar(int i, double a)
{
	if (i < 1 || i > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		for (int j = 1; j <= columns; j++) content[i - 1][j - 1] = content[i - 1][j - 1] * a;
	}
}

void Matrix::substract(int i, int j, double a)
{
	if (i < 1 || i > rows || j < 1 || j > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		for (int k = 1; k <= columns; k++)
			content[i - 1][k - 1] = content[i - 1][k - 1] - a * content[j - 1][k - 1];
	}
}

void Matrix::print()
{
	for (int i = 1; i <= rows; i++)
	{
		for (int j = 1; j < columns; j++) cout << this->getelement(i, j) << " ";
		cout << this->getelement(i, columns) << endl;
	}
}

Matrix operator* (const Matrix& A, const Matrix& B)
{
	if (A.columns != B.rows)
	{
		throw "\n GIVEN MATRICES ARE NOT MULTIPLIABLE.\n";
	}
	int m = A.rows;
	int n = A.columns;
	int p = B.columns;
	Matrix C(m, p);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < p; j++)
		{
			float z = 0;
			for (int k = 0; k < n; k++)
			{
				float x = A.content[i][k];
				float y = B.content[k][j];
				z += (x * y);
			}
			C.content[i][j] = z;
		}
	}
	return C;
}

Vector operator* (const Matrix& A, const Vector& X)
{
	if (A.columns != X.rows)
	{
		throw "\n GIVEN MATRICES ARE NOT MULTIPLIABLE.\n";
	}
	int m = A.rows;
	int n = A.columns;
	Vector Y(m);
	for (int i = 0; i < m; i++)
	{
		float b = 0;
		for (int k = 0; k < n; k++)
		{
			float a = A.content[i][k];
			float x = X.content[k][0];
			b += (a * x);
		}
		Y.content[i][0] = b;
	}
	return Y;
}

Vector operator-(const Vector& X, const Vector& Y0)
{
	Vector New = X;
	int N = New.getrows();
	Vector Y = Y0;
	for (int i = 0; i < N; i++)
	{
		New.content[i][0] = New.content[i][0] - Y.content[i][0];
	}
	return New;
}

double operator* (const Vector& X, const Vector& Y)
{
	double s = 0;
	for (int i = 0; i < X.rows; i++)
	{
		s += X.content[i][0] * Y.content[i][0];
	}
	return s;
}

Vector operator* (const double& a, const Vector& X)
{
	Vector Y(X.rows);
	for (int i = 0; i < X.rows; i++)
	{
		Y.content[i][0] = a * X.content[i][0];
	}
	return Y;
}

Vector operator+(const Vector& X, const Vector& Y0)
{
	Vector New = X;
	int N = New.getrows();
	Vector Y = Y0;
	for (int i = 0; i < N; i++)
	{
		New.content[i][0] += Y.content[i][0];
	}
	return New;
}

Matrix Matrix::T()
{
	Matrix T(columns, rows);
	for (int i = 1; i <= rows; i++)
	{
		for (int j = 1; j <= columns; j++)
		{
			T.setelement(j, i, this->content[i - 1][j - 1]);
		}
	}
	return T;
}