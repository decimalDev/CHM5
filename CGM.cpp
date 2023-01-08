#include "CGM.h"

Vector CGM(Lent& A, Vector& B, Vector& X0, double& quando, double eps)
{
	int N = A.get_N();
	Vector R = B - A.Prod(X0);
	Vector P = R;
	Vector X = X0;
	double a, b, rr;
	Vector Ap(N);
	Vector Buf(N);
	int iterationes = 0;
	do
	{
		Ap = A.Prod(P);
		rr = R * R;
		a = rr / (P * Ap);
		X = X + a * P;
		R = R - a * Ap;
		b = R * R / rr;
		P = R + b * P;
		iterationes++;
	}
	while (R.norm() >= eps);
	quando = iterationes;
	return X;
}

Vector PCGM(Lent& A, Vector& B, Vector& X0, double eps, int m, bool flag)
{
	int N = A.get_N();
	Vector R = B - A.Prod(X0);
	Vector S = mJacobi(A, R, X0, m);
	Vector P = S;
	Vector X = X0;
	double a, b, rs;
	Vector Ap(N);
	Vector Buf(N);
	int iterationes = 0;
	do
	{
		Ap = A.Prod(P);
		rs = R * S;
		a = rs / (P * Ap);
		X = X + a * P;
		R = R - a * Ap;
		S = mJacobi(A, R, X0, m);
		b = R * S / rs;
		P = S + b * P;
		iterationes++;
	} while (R.norm() >= eps);
	if(flag == false) cout << "PCGM with m = " << m << " requires " << iterationes << " iteration(s)." << endl;
	else
	{
		fstream fout;
		string FileName = "D://PCGM(m).txt";
		fout.open(FileName, ios::app);
		fout << endl << "___________________" << endl << m << "\t" << iterationes << endl;
		fout.close();
	}
	return X;
}

void CGM_as_func_of_q(int N, int l, double eps)
{
	cout << "CGM as func of (q,w)" << endl;
	//варьируем Q и w
	double quando;

	fstream fout;
	string FileName = "D://CGM(q).txt";
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
		CGM(A[q], B[q], X0, quando, eps);
		fout << Q[q] << "\t" << quando << endl;
		cout << Q[q] << "\t" << quando << endl;

	}
	fout.close();
}
