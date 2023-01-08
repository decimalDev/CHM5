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

#pragma once

using namespace std;

#pragma once

Vector CGM(Lent& A, Vector& B, Vector& X0, double& quando, double eps);

Vector PCGM(Lent& A, Vector& B, Vector& X0, double eps, int m, bool flag);

void CGM_as_func_of_q(int N, int l, double eps);
