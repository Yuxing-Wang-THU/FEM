#include <PhysicsEngine.h>
#include <Eigen\Dense>
#include <unsupported\Eigen\CXX11\Tensor>
#include<iostream>
#include <cmath>

// Mesh函数
Eigen::MatrixXi get_mesh(int element_num, int N_x, int N_y) {
    Eigen::MatrixXi element_idx(element_num, 3);
    element_idx.setZero();
    int cnt = 0;
    for (int j = 0; j < N_y - 1; j++) {
        for (int i = 0; i < N_x - 1; i++) {
            int idx = j * N_x + i;
            element_idx(cnt, 0) = idx;
            element_idx(cnt, 1) = idx + 1;
            element_idx(cnt, 2) = idx + N_x;
            cnt++;
            element_idx(cnt, 0) = idx + 1;
            element_idx(cnt, 1) = idx + N_x + 1;
            element_idx(cnt, 2) = idx + N_x;
            cnt++;
        }
    }
    return element_idx;
}

// 三角形面积
float get_current_area(Eigen::MatrixXd points) {
    float area = 0.5 * ((float)points(0, 0) * ((float)points(1, 1) - (float)points(2, 1)) + (float)points(1, 0) * (float)(points(2, 1) - (float)points(0, 1))
        + (float)points(2, 0) * ((float)points(0, 1) - (float)points(1, 1)));
    return area;
}

void set_zero_3dimtensor(Eigen::Tensor<Eigen::MatrixXd, 3>& tensor) {
    int shape = tensor.dimensions()[0];
    for (int i = 0; i < shape; i++) {
        tensor(i).setZero();
        }
}

// 计算D_m
void calculate_D_m(const Eigen::MatrixXd& node_pos, const Eigen::MatrixXi& element_idx, int element_num, Eigen::Tensor<Eigen::MatrixXd, 3>& D_m) {
    for (int ie = 0; ie < element_num; ie++) {
        Eigen::MatrixXd p0 = node_pos.row(element_idx(ie, 0));
        Eigen::MatrixXd p1 = node_pos.row(element_idx(ie, 1));
        Eigen::MatrixXd p2 = node_pos.row(element_idx(ie, 2));
        Eigen::MatrixXd dX(2, 2);
        dX << p1(0) - p0(0), p2(0) - p0(0),
              p1(1) - p0(1), p2(1) - p0(1);
        D_m(ie) = dX.inverse(); 
    }
}

// 计算dDsdx
// D_s = [[X2-X0  X4-X0
//        [X3-X1  X5-X1]]
// x1=[X0,X1],x2=[X2,X3],x3=[X4,X5]
void calculate_D_s(Eigen::Tensor<Eigen::MatrixXd, 3>& dDsdx) {
    Eigen::Matrix2d dX0, dX1, dX2, dX3, dX4, dX5;
    dX0 << -1, -1,
            0, 0;
    dDsdx(0) = dX0;
    dX1 << 0, 0,
          -1, -1;
    dDsdx(1) = dX1;
    dX2 << 1, 0,
           0, 0;
    dDsdx(2) = dX2;
    dX3 << 0, 0,
           1, 0;
    dDsdx(3) = dX3;
    dX4 << 0, 1,
           0, 0;
    dDsdx(4) = dX4;
    dX5 << 0, 0,
           0, 1;
    dDsdx(5) = dX5;
}

void implicit_euler_stvk(Eigen::MatrixXd& D_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, Eigen::Tensor<Eigen::MatrixXd, 3>& dDsdx,
Eigen::Tensor<Eigen::MatrixXd, 3>& dF, Eigen::Tensor<Eigen::MatrixXd, 3>& dE, Eigen::Tensor<Eigen::MatrixXd, 3>& dP, Eigen::Tensor<Eigen::MatrixXd, 3>& dH,
double init_area, double lam, double mu) {

    for (int i = 0; i < 6; i++) {
        dF(i) = dDsdx(i) * D_m;
        dE(i) = (dF(i).transpose() * F + F.transpose() * dF(i)) * 0.5;
        dP(i) = dF(i) * (2 * mu * strain + lam * strain.trace() * Eigen::MatrixXd::Identity(2, 2));
        dP(i) += F * (2 * mu * dE(i) + lam * dE(i).trace() * Eigen::MatrixXd::Identity(2, 2));
        dH(i) = -init_area * dP(i) * D_m.transpose();
    }
}

void implicit_euler_Neo(Eigen::MatrixXd& D_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, Eigen::Tensor<Eigen::MatrixXd, 3>& dDsdx,
    Eigen::Tensor<Eigen::MatrixXd, 3>& dF, Eigen::Tensor<Eigen::MatrixXd, 3>& dP, Eigen::Tensor<Eigen::MatrixXd, 3>& dH,
    double init_area, double lam, double mu) {
    
    Eigen::Matrix2d F_inv = F.inverse();

    double logJ = std::log(std::max(F.determinant(), 0.01));
    
    for (int i = 0; i < 6; i++) {
        dF(i) = dDsdx(i) * D_m;
        dP(i) = mu * dF(i) + (mu - lam * logJ) * F_inv.transpose() * dF(i).transpose() * F_inv.transpose();
        dP(i) += lam * (F_inv * dF(i)).trace() * F_inv.transpose();
        dH(i) = -init_area * dP(i) * D_m.transpose();
    }
}


Eigen::MatrixXd solver(int krow, int node_num, double mass, double dt, Eigen::MatrixXd& A, Eigen::MatrixXd& K, Eigen::VectorXd &b, Eigen::MatrixXd& node_vel, Eigen::MatrixXd& nodal_force) {
    // A
    for (int i = 0; i < krow; i++)
    {
        for (int j = 0; j < krow; j++)
        {
            if (i == j)
            {
                A(i, j) = mass - K(i, j) * dt * dt;
            }
            else
            {
                A(i, j) = -K(i, j) * dt * dt;
            }
        }
    }

    // b
    for (int k = 0; k < node_num; k++)
    {
        b(k * 2 + 0) = mass * node_vel(k, 0) + dt * nodal_force(k, 0);
        b(k * 2 + 1) = mass * node_vel(k, 1) + dt * nodal_force(k, 1);
    }

    Eigen::MatrixXd X = A.inverse() * b;
    return X;
}

// 2D形变函数
Eigen::MatrixXd deformation(Eigen::MatrixXd points, float theta, int tx, int ty, float sx, float sy) {
    double rad = M_PI * theta / 180.0;
    // 3x3旋转矩阵
    /*
     [[cos(rad),-sin(rad),0
       sin(rad),cos(rad),0
       0       ,0       ,1
     ]]
    */
    Eigen::Matrix3d rotation_matrix;
    rotation_matrix(0, 0) = cos(rad);
    rotation_matrix(0, 1) = -sin(rad);
    rotation_matrix(0, 2) = 0;
    rotation_matrix(1, 0) = sin(rad);
    rotation_matrix(1, 1) = cos(rad);
    rotation_matrix(1, 2) = 0;
    rotation_matrix(2, 0) = 0;
    rotation_matrix(2, 1) = 0;
    rotation_matrix(2, 2) = 1;

    // 平移矩阵
    Eigen::Matrix3d translation_matrix;
    translation_matrix << 1, 0, tx,
        0, 1, ty,
        0, 0, 1;
    // 缩放矩阵
    Eigen::Matrix3d scale_matrix;
    scale_matrix << sx, 0, 0,
        0, sy, 0,
        0, 0, 1;

    // 形变矩阵
    Eigen::Matrix3d transform_matrix = translation_matrix * rotation_matrix * scale_matrix;

    // n x 3大小的坐标矩阵
    Eigen::MatrixXd rts(points.rows(), 3);
    //block(i,j,p,q)，提取块大小为(p,q),起始于(i,j)
    rts.block(0, 0, points.rows(), 2) = points;
    //最后一列搞成1
    rts.col(2).setOnes();
    //返回形变后点的坐标
    Eigen::MatrixXd result = (rts * transform_matrix.transpose()).leftCols(2);
    return result;
}

// Linear本构模型
Eigen::Matrix2d LinearFemStep(Eigen::MatrixXd& node_deformed, Eigen::Matrix2d B_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, double init_area, double lam, double mu) {
    // Step 1. Calculate D_s (2 x 2)
    Eigen::Matrix2d D_s;
    D_s(0, 0) = node_deformed(1, 0) - node_deformed(0, 0);
    D_s(0, 1) = node_deformed(2, 0) - node_deformed(0, 0);
    D_s(1, 0) = node_deformed(1, 1) - node_deformed(0, 1);
    D_s(1, 1) = node_deformed(2, 1) - node_deformed(0, 1);

    // Step 3. Calculate Deformation gradient
    F = D_s * B_m;

    // Step 4. Calculate PK1 stress of Linear Elasticticy
    strain = (F + F.transpose()) * 0.5 - Eigen::Matrix2d::Identity();
    double doubleInner = strain(0, 0) * strain(0, 0) + strain(1, 0) * strain(1, 0) + strain(0, 1) * strain(0, 1) + strain(1, 1) * strain(1, 1);
    double energy = doubleInner * mu + lam * 0.5 * pow(strain.trace(), 2);

    Eigen::Matrix2d piola = mu * (F + F.transpose() - 2 * Eigen::Matrix2d::Identity()) + lam * (F(0, 0) - 1 + F(1, 1) - 1) * Eigen::Matrix2d::Identity();

    // Step 5. Calculate Force
    Eigen::Matrix2d H = -init_area * piola * B_m.transpose();

    return H;
}

// STVK本构模型
Eigen::Matrix2d StvkFemStep(Eigen::MatrixXd& node_deformed, Eigen::Matrix2d B_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, double init_area, double lam, double mu) {
    // Step 1. Calculate D_s (2 x 2)
    Eigen::Matrix2d D_s;
    D_s(0, 0) = node_deformed(1, 0) - node_deformed(0, 0);
    D_s(0, 1) = node_deformed(2, 0) - node_deformed(0, 0);
    D_s(1, 0) = node_deformed(1, 1) - node_deformed(0, 1);
    D_s(1, 1) = node_deformed(2, 1) - node_deformed(0, 1);

    // Step 3. Calculate Deformation gradient
    F = D_s * B_m;

    // Step 4. Calculate PK1 stress 
    strain = ((F.transpose() * F) - Eigen::Matrix2d::Identity()) * 0.5;
    double doubleInner = strain(0, 0) * strain(0, 0) + strain(1, 0) * strain(1, 0) + strain(0, 1) * strain(0, 1) + strain(1, 1) * strain(1, 1);
    double energy = doubleInner * mu + lam * 0.5 * pow(strain.trace(), 2);

    Eigen::Matrix2d piola = F * (2 * mu * strain + lam * strain.trace() * Eigen::Matrix2d::Identity());

    // Step 5. Calculate Force
    Eigen::Matrix2d H = -init_area * piola * B_m.transpose();

    return H;
}

// CoRotated本构模型
Eigen::Matrix2d CoRotatedFemStep(Eigen::MatrixXd& node_deformed, Eigen::Matrix2d B_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, double init_area, double lam, double mu) {

    // Step 1. Calculate D_s (2 x 2)
    Eigen::Matrix2d D_s;
    D_s(0, 0) = node_deformed(1, 0) - node_deformed(0, 0);
    D_s(0, 1) = node_deformed(2, 0) - node_deformed(0, 0);
    D_s(1, 0) = node_deformed(1, 1) - node_deformed(0, 1);
    D_s(1, 1) = node_deformed(2, 1) - node_deformed(0, 1);

    // Step 3. Calculate Deformation gradient
    F = D_s * B_m;

    Eigen::JacobiSVD<Eigen::Matrix2d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix2d U = svd.matrixU();
    Eigen::Vector2d sigma = svd.singularValues();
    Eigen::Matrix2d Vt = svd.matrixV().transpose();

    // Get polar decomposition
    Eigen::Matrix2d R = U * Vt;
    Eigen::Matrix2d S = R.inverse() * F;
    strain = S - Eigen::Matrix2d::Identity();
    // Get PK1 stress
    // piola = 2µ(F − R) + λtr(RTF − I)R
    Eigen::Matrix2d RTF_I = R.transpose() * F - Eigen::Matrix2d::Identity();
    Eigen::Matrix2d piola = 2 * mu * (F - R) + lam * RTF_I.trace() * R;

    // Step 5. Calculate Force
    Eigen::Matrix2d H = -init_area * piola * B_m.transpose();
    return H;
}

// Neohookean本构模型
Eigen::Matrix2d NeohookeanFemStep(Eigen::MatrixXd& node_deformed, Eigen::Matrix2d B_m, Eigen::Matrix2d& F, Eigen::Matrix2d& strain, double init_area, double lam, double mu) {

    // Step 1. Calculate D_s (2 x 2)
    Eigen::Matrix2d D_s;
    D_s(0, 0) = node_deformed(1, 0) - node_deformed(0, 0);
    D_s(0, 1) = node_deformed(2, 0) - node_deformed(0, 0);
    D_s(1, 0) = node_deformed(1, 1) - node_deformed(0, 1);
    D_s(1, 1) = node_deformed(2, 1) - node_deformed(0, 1);

    // Step 3. Calculate Deformation gradient
    F = D_s * B_m;

    // Calculate strain tensor of Linear Elasticity
     // Green strain tensor e = 1/2(F^TF-I)
    strain = (F.transpose() * F - Eigen::Matrix2d::Identity()) * 0.5;

    // Get Log J
    double logJ = log(fmax(F.determinant(), 0.01));
    // double log_I3 = std::log(std::max(F.determinant() * F.determinant(), 0.01));

    double I_1 = (F.transpose() * F).trace();
    // Calculate Energy density function for updating position

    // Ψ(I1, J) = µ/2(I1 - 2) - µlog(J) + λ/2log^2(J)
    double energy = 0.5 * mu * (I_1 - 2) - mu * logJ + 0.5 * lam * logJ * logJ;

    // Ψ(I1, J) = µ/2(I1 - 2) - µlog(J) + λ/2log^2(J)
    // energy = 0.5 * mu * (I_1 - log_I3 - 2) + lam/8 * log_I3 * log_I3;

    // Get PK1 stress
    // piola = µ(F - µF^-T) + λ log(J)F^-T
    Eigen::Matrix2d F_inv = F.inverse();
    Eigen::Matrix2d F_inv_T = F_inv.transpose();

    Eigen::Matrix2d piola = mu * (F - F_inv_T) + lam * logJ * F_inv_T;
    // piola = mu * (F - F_inv_T) + lam * 0.5 * log_I3 * F_inv_T;


    // Step 5. Calculate Force
    Eigen::Matrix2d H = -init_area * piola * B_m.transpose();
    return H;
}
