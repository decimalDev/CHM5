#include "Task I.h"

vector<double> NonDiagonal(int N, int l)
{
	vector<double> A;
	int card = N * (2 * l + 1) - l * (l + 1); //общее число элементов ленты
	A.resize(card);
	for (int i = N; i < card; i++)
	{
		int num = -1000 + rand() % (2000 + 1);
		double a = (double)num / 1000.0;
		A[i] = a;
	}
	return A;
}

int index(int i, int j, int N, int l)
{
	int d = i - j;
	int D = abs(d);
	int index;
	if (d == 0)
	{
		index = i - 1;
	}
	if (d < 0)
	{
		index = D * (2 * N - D) + j - D - 1;
	}
	if (d > 0)
	{
		index = N + (D - 1) * (2 * N - D) + i - D - 1;
	}
	return index;
}

double Lent::norm()
{
	double n = 0;
	for (int i = N; i < L.size(); i++)
	{
		if (abs(L[i]) > n) n = abs(L[i]);
	}
	return n;
}

bool Lent::IfDiagDom()
{
	vector<double> S = this->RowSum();
	for (int i = 0; i < N; i++)
	{
		if (L[i] <= S[i]) return false;
	}
	return true;
}

Lent::Lent(int N, int l)
{
	vector<double> A;
	double card = N * (2 * l + 1) - l * (l + 1); //общее число элементов ленты
	A.resize(card, 0.0);
	this->N = N;
	this->l = l;
	this->L = A;
}

Lent::Lent(int N, int l, vector<double>& L)
{
	this->N = N;
	this->l = l;
	this->L = L;
}

int Lent::get_N()
{
	return N;
}

int Lent::get_l()
{
	return l;
}

double Lent::get_element(int i, int j)
{
	return L[index(i, j, N, l)];
}

vector<double> Lent::RowSum()
{
	vector<double> S;
	for (int row = 1; row <= N; row++)
	{
		double s = 0;
		int left = max(1, row - l);
		int column = left;
		while(column < row)
		{
			s += abs(L[index(row, column, N, l)]);
			column++;
		}
		int right = min(N, row + l);
		column = row + 1;
		while(column <= right)
		{
			s += abs(L[index(row, column, N, l)]);
			column++;
		}
		S.push_back(s);
	}
	return S;
}

void Lent::print()
{
	for (int i = 1; i <= N; i++)
	{
		int left = max(1, i - l);
		int j = 1;
		while (j < left)
		{
			cout << "0.000 ";
			j++;
		}
		while (j <= i)
		{
			cout << L[index(i, j, N, l)] << " ";
			j++;
		}
		int right = min(N, i + l);
		while(j <= right)
		{
			cout << L[index(i, j, N, l)] << " ";
			j++;
		}
		while (j <= N)
		{
			cout << "0.000 ";
			j++;
		}
		cout << endl;
	}
}

Vector Lent::Prod(Vector& X)
{
	Vector Y(N);
	for (int i = 1; i <= N; i++)
	{
		double s = 0;
		int left = max(1, i - l);
		int right = min(N, i + l);
		for (int j = left; j <= right; j++) { s += L[index(i, j, N, l)] * X.getelement(j); }
		Y.setelement(i, s);
	}
	return Y;
}

Lent::Lent(const Lent& copia)
{
	N = copia.N;
	l = copia.l;
	L = copia.L;
}

Lent& Lent::operator=(const Lent& copia)
{
	N = copia.N;
	l = copia.l;
	L = copia.L;
	return *this;
}

Lent Lent::T()
{
	Lent A(N, l);
	for (int i = 1; i <= N; i++)
	{
		int left = max(1, i - l);
		int right = min(N, i + l);
		int j = left;
		while (j <= right)
		{
			A.L[index(j, i, N, l)] = this->L[index(i, j, N, l)];
			j++;
		}
	}
	return A;
}

Lent operator*(const Lent& A, const Lent& B)
{
	int N = A.N;
	int l = A.l;
	Lent P(N, 2*l);
	for (int i = 1; i <= N; i++)
	{
		int left_i = max(1, i - l);
		int right_i = min(N, i + l);
		for (int j = max(1, i - 2*l); j <= min(N, 2*l + i); j++)
		{
				int left_j = max(1, j - l);
				int right_j = min(N, j + l);
				int left = max(left_i, left_j);
				int right = min(right_i, right_j);
				int k = left;
				double s = 0;
				while (k <= right)
				{
					s += A.L[index(i, k, N, l)] * B.L[index(k, j, N, l)];
					k++;
				}
				P.L[index(i, j, N, l)] = s;
		}
	}
	return P;
}


vector<Vector> Generate_Lents(vector<Lent>& A, Vector& X, vector<double> Q, int N, int l) //генерируем Q.size() штук l-диагональных положительно определённых симметричных квадратных матриц порядка N
{
	A = {};
	int card = Q.size();
	vector<double> v = NonDiagonal(N, l);  //недиагональные элементы ленты
	cout << endl;
	Lent U(N, l, v); //строим исходную матрицу
	vector<double> S = U.RowSum(); //считаем суммы модулей недиагональных элементов в каждой строке
	vector<vector<double>> vv;
	for (int i = 0; i < card; i++)
	{
		vv.push_back(v);
	}
	vector<Lent> A0; //A0 -- исходные ленточные матрицы, A -- после приведения к оптимальному виду
	vector<Vector> B0; //исходный вектор значений и вектор значений после домножения на A^T
	vector<Vector> B;
	for (int i = 0; i < card; i++)
	{
		for (int j = 0; j < N; j++)
		{
			vv[i][j] = Q[i] * S[j]; //диагональные элементы
		}
		Lent Aa(N, l, vv[i]); //матрица A_i из заданий
		A0.push_back(Aa);
		cout << "\tA(" << i + 1 << "): " << endl;
		//A0[i].print();  //сходится: недиагональные элементы матриц A1-A3 совпадают, диагональные пропорциональны с требуемой пропорциональностью
		Vector b = A0[i].Prod(X); //вектор значений получаем как произведение Ax
		B0.push_back(b);
		cout << "\tB(" << i + 1 << "): " << endl;
		// (int i = 1; i <= M; i++) cout << b.getelement(i) << " ";
		cout << endl;
		Lent At(N, l); //транспонированные матрицы
		At = A0[i].T();
		Lent a(N, 2 * l); //матрица A^T * A
		a = At * A0[i];
		A.push_back(a);
		cout << "\tA*(" << i + 1 << "): " << endl;
		//A[i].print();
		if (A[i].IfDiagDom() == false) cout << "NOT DIAGONALLY DOMINANT!" << endl; //проверяем диагональное преобладание -- во многих случаях нарушается
		B.push_back(At.Prod(B0[i])); //умножаем исходное B на A^T
		//cout << "\tB*(" << i + 1 << "): " << endl;
		//for (int j = 1; j <= M; j++) cout << B[i].getelement(j) << " ";
		cout << endl;
		cout << endl;
	}
	return B;
}