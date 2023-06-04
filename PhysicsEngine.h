#define _USE_MATH_DEFINES 
#ifndef _DEFORMATION_H
#define _DEFORMATION_H

#include <unsupported\Eigen\CXX11\Tensor>
#include <Eigen\Dense>
#include <iostream>
#include <cmath>

Eigen::MatrixXi get_mesh(int element_num, int N_x, int N_y);

float get_current_area(Eigen::MatrixXd points);

void calculate_D_s(Eigen::Tensor<Eigen::MatrixXd, 3>& dDsdx);

void set_zero_3dimtensor(Eigen::Tensor<Eigen::MatrixXd, 3>& tensor);

void calculate_D_m(const Eigen::MatrixXd& node_pos, const Eigen::MatrixXi& element_idx, int element_num, Eigen::Tensor<Eigen::MatrixXd, 3>& D_m);

void implicit_euler_stvk(Eigen::MatrixXd& D_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, Eigen::Tensor<Eigen::MatrixXd, 3>& dDsdx,
	Eigen::Tensor<Eigen::MatrixXd, 3>& dF, Eigen::Tensor<Eigen::MatrixXd, 3>& dE, Eigen::Tensor<Eigen::MatrixXd, 3>& dP, Eigen::Tensor<Eigen::MatrixXd, 3>& dH,
	double init_area, double lam, double mu);

void implicit_euler_Neo(Eigen::MatrixXd& D_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, Eigen::Tensor<Eigen::MatrixXd, 3>& dDsdx,
	Eigen::Tensor<Eigen::MatrixXd, 3>& dF, Eigen::Tensor<Eigen::MatrixXd, 3>& dP, Eigen::Tensor<Eigen::MatrixXd, 3>& dH,
	double init_area, double lam, double mu);

Eigen::MatrixXd solver(int krow, int node_num, double mass, double dt, Eigen::MatrixXd& A, Eigen::MatrixXd& K, Eigen::VectorXd& b, Eigen::MatrixXd& node_vel, Eigen::MatrixXd& nodal_force);

Eigen::MatrixXd deformation(Eigen::MatrixXd points, float theta, int tx, int ty, float sx, float sy);

Eigen::Matrix2d LinearFemStep(Eigen::MatrixXd& node_deformed, Eigen::Matrix2d B_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, double init_area, double lam, double mu);

Eigen::Matrix2d StvkFemStep(Eigen::MatrixXd& node_deformed, Eigen::Matrix2d B_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, double init_area, double lam, double mu);

Eigen::Matrix2d CoRotatedFemStep(Eigen::MatrixXd& node_deformed, Eigen::Matrix2d B_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, double init_area, double lam, double mu);

Eigen::Matrix2d NeohookeanFemStep(Eigen::MatrixXd& node_deformed, Eigen::Matrix2d B_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, double init_area, double lam, double mu);

#endif