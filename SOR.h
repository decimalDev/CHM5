#include<iostream>
#include<fstream>
#include<vector>
#include<stdio.h>
#include<cmath>
#include <iomanip>
#include <stdlib.h>
#include "Matrix algebra.h"
#include "Task I.h"

#pragma once

using namespace std;

Vector SOR(Lent& A, Vector& B, Vector& X0, double eps, double& quando, double omega);

void SOR_as_func_of_q_w(int N, int l, double eps);
