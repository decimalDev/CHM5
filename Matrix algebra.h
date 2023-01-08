#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

#pragma once

class Vector;

// ����� ������� ��������������� ����� ������: ������ ����� (rows), ������ �������� (columns) � ��������� �������� ��������� (content)
class Matrix
{
protected:
	int rows;
	int columns;
	vector<vector<double>> content;

public:
	int getrows();

	int getcolumns();

	void setcontent(int m, int n); // �� ������� content -- ������� ������� ���������������� �������

	//�����, ����������� ��������� ��������� ������ ������ ������
	void setelement(int i, int j, double a);

	//�����, ����������� �������� �� ������ ������
	double getelement(int i, int j);

	//�����, �������������� ������������ �����
	void permutation(int i, int j);

	//�����, �������������� ��������� ����� �� ������
	void scalar(int i, double a);

	//�����, ���������� ����������� �� a j-� ������ �� i-��
	void substract(int i, int j, double a);

	//����������������
	Matrix T();

	//��������� ����������� ������ ������� 3x3
	Matrix();

	Matrix(int m, int n);

	Matrix(int m, int n, double k); //������ ������� m x n � ���������� k �� ������� ���������

	Matrix(vector<vector<double>>& C);

	void print();

	friend Matrix operator* (const Matrix& A, const Matrix& B);

	friend Vector operator* (const Matrix& A, const Vector& X);

};

class Vector : public Matrix
{
public:
	void setelement(int i, double a);

	double getelement(int i);

	double norm();

	Vector();

	Vector(int m);
	Vector(const Vector& copia);

	Vector(vector<double>& C);

	//friend Vector operator* (const RotationMatrix& A, const Vector& X);

	Vector& operator= (const Vector& v);

	friend Vector operator* (const Matrix& A, const Vector& X);
	friend Vector operator- (const Vector& X, const Vector& Y);
	friend double operator* (const Vector& X, const Vector& Y);
	friend Vector operator* (const double& a, const Vector& X);
	friend Vector operator+ (const Vector& X, const Vector& Y);

};

