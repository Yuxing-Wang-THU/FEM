#ifndef _MAIN_H
#define _MAIN_H
#include <iostream>
#include <set>

using namespace std;

// Visualization
int RESOLUTION_X = 1280;
int RESOLUTION_Y = 480;

// Terrain
int L_U_X = 0;
int L_U_Y = 0;
int R_D_X = RESOLUTION_X;
int R_D_Y = 30;

// Iteration
int ITER = 0;
int MAX_ITER = 500;

// Young's modulus
double Y = 10000;

// Poisson's ratio
double V = 0.1;

// Lame
double MU = Y / (2 * (1 + V));
double LAM = Y * V / (1 + V) / (1 - 2 * V);

// Time horizon
double DT = 0.0001;

// Mass
double POINT_MASS = 1.0;

// Cube
int VOXEL_X = 8;
int VOXEL_Y = 3;
int CUBE_LEN = 1;
int INIT_X = 0;
int INIT_Y = 4;

int N_X = VOXEL_X + 1;
int N_Y = VOXEL_Y + 1;
int NODE_NUM = N_X * N_Y;
int ELEMENT_NUM = 2 * VOXEL_X * VOXEL_Y;

set<int> EMPTY_INDEX {2, 3, 6, 7, 10, 11, 14, 15};
//set<int> EMPTY_INDEX {};

string MODEL = "Neo";
string TIME_INTE = "IMP";

#endif