#include<iostream>
#include <fstream>
#include<vector>
#include<stdio.h>
#include<cmath>
#include <iomanip>
#include <stdlib.h>
#include "Matrix algebra.h"
#include "Task I.h"

#pragma once

using namespace std;

Vector Jacobi(Lent& A, Vector& B, Vector& X, double& quando, double eps);

Vector mJacobi(Lent& A, Vector& B, Vector& X0, int m);

void Jacobi_as_func_of_q(int N, int l, double eps);
