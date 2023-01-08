#include "Jacobi.h"

Vector Jacobi(Lent& A, Vector& B, Vector& X0, double& quando, double eps)
{
	int N = X0.getrows();
	int l = A.get_l();

	Vector X = X0;
	Vector Xnew(N);
	Vector Buf(N);
	Vector r(N);
	int iterationes = 0;
	do
	{
		for (int i = 1; i <= N; i++)
		{
			double x = 0;
			int left = max(1, i - l);
			int right = min(i + l, N);
			int j = left;
			while(j < i)
			{
				x += A.get_element(i, j) * X.getelement(j);
				j++;
			}
			j = i + 1;
			while (j <= min(i + l, N))
			{
				x += A.get_element(i, j) * X.getelement(j);
				j++;
			}
			x = (B.getelement(i) - x) / A.get_element(i, i);
			Xnew.setelement(i, x);
		}
		Buf = X;
		X = Xnew;
		r = A.Prod(X) - B;
		iterationes++;
	}
	while (r.norm() >= eps);
	quando = iterationes;
	return X;
}

Vector mJacobi(Lent& A, Vector& B, Vector& X0, int m)
{
	int N = X0.getrows();
	int l = A.get_l();

	Vector X = X0;
	Vector Xnew(N);
	int p = 0;
	do
	{
		for (int i = 1; i <= N; i++)
		{
			double x = 0;
			int left = max(1, i - l);
			int right = min(i + l, N);
			int j = left;
			while (j < i)
			{
				x += A.get_element(i, j) * X.getelement(j);
				j++;
			}
			j = i + 1;
			while (j <= min(i + l, N))
			{
				x += A.get_element(i, j) * X.getelement(j);
				j++;
			}
			x = (B.getelement(i) - x) / A.get_element(i, i);
			Xnew.setelement(i, x);
		}
		X = Xnew;
		p++;
	} while (p <= m);
	return X;
}

void Jacobi_as_func_of_q(int N, int l, double eps)
{
	cout << "Jacobi as func of q" << endl;
	//варьируем Q и w
	double quando;

	fstream fout;
	string FileName = "D://Jacobi(q).txt";
	fout.open(FileName, ios::out);

	Vector X(N); //генерируем вектор X
	cout << "\tX: " << endl;
	for (int i = 1; i <= N; i++)
	{
		int num = -1000 + rand() % (2000 + 1);
		double a = (double)num / 1000.0;
		X.setelement(i, a);
	}
	cout << endl;

	Vector X0(N); //начальное приближение -- нулевой вектор

	//генерируем A, Q и B
	vector<Lent> A;
	vector<double> Q;
	double i = 1.5;
	while (i < 10)
	{
		Q.push_back(i);
		i += 0.5;
	}
	vector<Vector> B = Generate_Lents(A, X, Q, N, l); //матрицы
	for (int q = 0; q < Q.size(); q++)
	{
			Jacobi(A[q], B[q], X0, quando, eps);
			fout << Q[q] << "\t" << quando << endl;
			cout << Q[q] << "\t" << quando << endl;

	}
	fout.close();
}