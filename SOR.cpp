#include "SOR.h"

Vector SOR(Lent& A, Vector& B, Vector& X0, double eps, double& quando, double w)
{
	int N = A.get_N();
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
			while (j < i)
			{
				x += A.get_element(i, j) * Xnew.getelement(j);
				j++;
			}
			j = i + 1;
			while (j <= min(i + l, N))
			{
				x += A.get_element(i, j) * X.getelement(j);
				j++;
			}
			x = (B.getelement(i) - x) * w / A.get_element(i, i);
			x += (1 - w) * X.getelement(i);
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

void SOR_as_func_of_q_w(int N, int l, double eps)
{
	cout << "SOR as func of (q,w)" << endl;
	//варьируем Q и w
	double quando;

	fstream fout;
	string FileName = "D://SOR(q,w).txt";
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
		double o = 0.2;
		while (o < 2)
		{
			cout << Q[q] << "\t" << o << endl;
			SOR(A[q], B[q], X0, eps, quando, o);
			fout << Q[q] << "\t" << o << "\t" << quando << endl;
			o += 0.1;
		}
	}
	fout.close();
}