#include<iostream>
#include <fstream>
#include<vector>
#include<stdio.h>
#include<cmath>
#include <iomanip>
#include <stdlib.h>
#include "Matrix algebra.h"
#include "Task I.h"
#include "Jacobi.h"
#include "SOR.h"
#include "CGM.h"
#include "Matrix algebra.cpp"
#include "Task I.cpp"
#include "Jacobi.cpp"
#include "SOR.cpp"
#include "CGM.cpp"

#pragma once

int M = 1000; //размерность матриц
int k = 40; //ширина ленты
double eps = 0.000001;
double w = 1.0; //параметр omega для метода SOR
double q1 = 1.1;
double q2 = 2.0;
double q3 = 10.0;
double m = 10; //число шагов в предобуславливании в PCGM

using namespace std;

int main1()
{
	cout.width(14);
	cout.setf(ios::left);
	cout.precision(10);
	cout << endl;
	srand(time(NULL));

	//Jacobi_as_func_of_q(M, k, eps);

	//SOR_as_func_of_q_w(M, k, eps); //строим график зависимости числа итераций метода SOR от q и w

	//CGM_as_func_of_q(M, k, eps);
	

	//vector<double> testv = { 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1, 10.1, 11.1, 12.1, 13.1, 14.1, 15.1, 16.1, 17.1, 18.1, 19.1 }; //проверим работу основных функций на простой матрице
	//Lent testLent(5, 2, testv);

	//cout << "Test matrix: " << endl;
	//testLent.print(); //вывод верный
	//cout << endl;
	//Lent testLentT = (testLent.T());
	///cout << "Testing transpose operation: " << endl;
	//testLentT.print(); //вывод верный
	//cout << endl;

	//Lent testP = testLentT * testLent;
	//cout << "Testing A^T * A product: " << endl;
	//testP.print(); //вывод верный
	//cout << endl;

	//cout << "Testing RowSum:" << endl;
	//vector<double> testRS = testP.RowSum();
	//for (double x : testRS) cout << x << " "; //вывод верный
	//cout << endl;
	// 
	//task1
	cout << "task1"<<endl;
	Vector X(M); //генерируем вектор X
	cout << "\tX: " << endl;
	for (int i = 1; i <= M; i++)
	{
		int num = -1000 + rand() % (2000 + 1);
		double a = (double)num / 1000.0;
		X.setelement(i, a);
		//cout << X.getelement(i) << " ";
	}
	cout << endl;

	//генерируем матрицы A1, A2, A3 и вектор B
	vector<Lent> A;
	vector<double> Q = { q1, q2, q3 };
	vector<Vector> B = Generate_Lents(A, X, Q, M, k); //три базовые матрицы

	//task2
	cout << "task2"<<endl;
	//МЕТОД ЯКОБИ
	Vector X0(M); //в качестве начального приближения используем нулевой вектор
	vector<Vector> XJacobi; //вектор решений
	cout << endl;
	for (int i = 0; i < 3; i++)
	{
		double quando;
		XJacobi.push_back(Jacobi(A[i], B[i], X0, quando, eps));
		cout << "X_{Jacobi}(" << i + 1 << "): " << endl;
		//for (int j = 1; j <= M; j++)
		//{
		//	cout << XJacobi[i].getelement(j) << " ";
		//}
		if ((X - XJacobi[i]).norm() < eps) cout << "OK"; //проверяем точность
		else cout << "!" << setw(7) << setprecision(4) << (X - XJacobi[i]).norm();
		cout << endl;
	}
	cout << endl;

	cout << "task3" << endl;
	//task3
	//SOR
	double quando;
	vector<Vector> XSOR; //вектор решений
	for (int i = 0; i < 3; i++)
	{
		XSOR.push_back(SOR(A[i], B[i], X0, eps, quando, w));
		cout << "X_{SOR}(" << i + 1 << "): " << endl;
		/*for (int j = 1; j <= M; j++)
		{
			cout << XSOR[i].getelement(j) << " ";
		}*/
		if ((X - XSOR[i]).norm() < eps) cout << "OK"; //проверяем точность
		else cout << "!" << setw(7) << setprecision(4) << (X - XSOR[i]).norm();
		cout << endl;
	}
	cout << endl;

	cout << "task4" << endl;
	//task4
	//CGM

	//vector<double> testCGM = { 4.0, 3.0, 1.0, 1.0 }; //тестируем CGM на примере из Википедии
	//Lent AtestCGM(2, 1, testCGM);
	//vector<double> ytestCGM = { 1.0, 2.0 };
	//Vector YtestCGM(ytestCGM);
	//vector<double> x0testCGM = { 2.0, 1.0 };
	//Vector X0testCGM(x0testCGM);
	//Vector XtestCGM = CGM(AtestCGM, YtestCGM, X0testCGM, 0.0001);
	//cout << "\t XtestCGM: " << endl; //работает!
	//cout << XtestCGM.getelement(1) << " " << XtestCGM.getelement(2) << endl;

	vector<Vector> XCGM; //вектор решений
	for (int i = 0; i < 3; i++)
	{
		XCGM.push_back(CGM(A[i], B[i], X0, quando, eps));
		cout << "X_{CGM}(" << i + 1 << "): " << endl;
		//for (int j = 1; j <= M; j++)
		//{
		//	cout << XCGM[i].getelement(j) << " ";
		//}
		if ((X - XCGM[i]).norm() < eps) cout << "OK"; //проверяем точность
		else cout << "!" << setw(7) << setprecision(4) << (X - XCGM[i]).norm();
		cout << endl;
	}
	cout << endl;

	vector<Vector> XPCGM; //вектор решений
	for (int i = 0; i < 3; i++)
	{
		XPCGM.push_back(PCGM(A[i], B[i], X0, eps, m, false));
		cout << "X_{PCGM}(" << i + 1 << "): " << endl;
		//for (int j = 1; j <= M; j++)
		//{
		//	cout << XPCGM[i].getelement(j) << " ";
		//}
		if ((X - XPCGM[i]).norm() < eps) cout << "OK"; //проверяем точность
		else cout << "!" << setw(7) << setprecision(4) << (X - XPCGM[i]).norm();
		cout << endl;
	}

	//варьируем число шагов в методе PCGM и смотрим, как меняется число итераций
	//for (int mu = 1; mu <= 30; mu++)
	//{
	//	PCGM(A[0], B[0], X0, eps, mu, true);
	//}

	return 0;
}