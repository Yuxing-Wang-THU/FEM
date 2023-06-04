#include <main.h>
#include <PhysicsEngine.h>
#include <matplotlibcpp.h>
#include <Eigen\Dense>
#include <unsupported\Eigen\CXX11\Tensor>
#include <vector>
#include <iostream>

namespace plt = matplotlibcpp;
using namespace std;

int main()
{   // 初始化PLT
    plt::figure_size(1280, 720);

    // Mesh
    Eigen::MatrixXi element_idx = get_mesh(ELEMENT_NUM, N_X, N_Y);
    
    // 位置，速度，节点力
    Eigen::MatrixXd node_pos = Eigen::MatrixXd::Zero(NODE_NUM, 2);
    Eigen::MatrixXd node_vel = Eigen::MatrixXd::Zero(NODE_NUM, 2);
    Eigen::MatrixXd nodal_force = Eigen::MatrixXd::Zero(NODE_NUM, 2);

    // 重力
    Eigen::MatrixXd gravity(1, 2);
    gravity << 0.0, -POINT_MASS * 9.8;

    // 设置初始位置
    for (int j = 0; j < N_Y; j++) {
        for (int i = 0; i < N_X; i++) {
            node_pos.row(j * N_X + i) << INIT_X + i * CUBE_LEN, INIT_Y + j * CUBE_LEN;
        }
    }
    
    // Element的初始面积
    Eigen::MatrixXd init_tri(3, 2);
    for (int i = 0; i < 3; i++) {
        init_tri.row(i) = node_pos.row(element_idx(0, i));
    }
    float init_area = get_current_area(init_tri);

    // D_m矩阵 
    Eigen::Tensor<Eigen::MatrixXd, 3> D_m(ELEMENT_NUM, 2, 2);
    calculate_D_m(node_pos, element_idx, ELEMENT_NUM, D_m);

    // 形变
    Eigen::MatrixXd position = deformation(node_pos, 0.0, 0, 0, (float)0.9, (float)1.0);

    // 调PLT
    std::map<std::string, std::string> kwargs;
    {
        for (int iele = 0; iele < ELEMENT_NUM; iele++) {
            if (EMPTY_INDEX.count(iele) > 0) {
                continue;
            }
            std::vector<double> element_x{}, element_y{};
            for (int point = 0; point < 3; point++) {
                element_x.push_back(position(element_idx(iele, point), 0));
                element_y.push_back(position(element_idx(iele, point), 1));
            }
            plt::fill(element_x, element_y, kwargs);
        }
        plt::draw();
        plt::xlim(0, 30);
        plt::ylim(0, 15);
        plt::tight_layout();
    }

    // 隐式方法
    // 构建线性方程组AX=b
    const int krow = NODE_NUM * 2;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(krow, krow);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(krow, krow);
    Eigen::VectorXd b(krow);
    b.setZero();

    // ds, dF, dP, dH矩阵 
    Eigen::Tensor<Eigen::MatrixXd, 3> dDsdx(6, 2, 2);
    Eigen::Tensor<Eigen::MatrixXd, 3> dF(6, 2, 2);
    Eigen::Tensor<Eigen::MatrixXd, 3> dP(6, 2, 2);
    Eigen::Tensor<Eigen::MatrixXd, 3> dH(6, 2, 2);
    Eigen::Tensor<Eigen::MatrixXd, 3> dE(6, 2, 2);

    //FEM_MAIN_LOOP
    while (ITER < MAX_ITER) {
        ITER += 1;
        
        // 清屏
        plt::clf();
        
        //调PLT
        {
            for (int iel = 0; iel < ELEMENT_NUM; iel++) {
                if (EMPTY_INDEX.count(iel) > 0) {
                    continue;
                }
                std::vector<double> element_x{}, element_y{};
                for (int point = 0; point < 3; point++) {
                    element_x.push_back(position(element_idx(iel, point), 0));
                    element_y.push_back(position(element_idx(iel, point), 1));
                }
                plt::fill(element_x, element_y, kwargs);
            }
            plt::draw();
            plt::xlim(0, 30);
            plt::ylim(0, 15);
            plt::tight_layout();
            plt::pause(0.01);
        }
        
        // Substep
        Eigen::Matrix2d F, strain;
        for (int i = 0; i < 30; i++) {
            nodal_force.setZero();
            for (int ie = 0; ie < ELEMENT_NUM; ie++) {
                if (EMPTY_INDEX.count(ie) > 0) {
                    continue;
                }
                // 计算节点力
                Eigen::MatrixXd element_pos(3, 2);
                for (int j = 0; j < 3; j++) {
                    element_pos.row(j) = position.row(element_idx(ie, j));
                }

                if (MODEL == "Linear") {
                    Eigen::Matrix2d H = LinearFemStep(element_pos, D_m(ie), F,strain, init_area, LAM, MU);
                    nodal_force.row(element_idx(ie, 1)) += H.col(0);
                    nodal_force.row(element_idx(ie, 2)) += H.col(1);
                    nodal_force.row(element_idx(ie, 0)) += (-H.col(0) - H.col(1));
                }

                else if (MODEL == "STVK") {
                    Eigen::Matrix2d H = StvkFemStep(element_pos, D_m(ie), F, strain, init_area, LAM, MU);
                    nodal_force.row(element_idx(ie, 1)) += H.col(0);
                    nodal_force.row(element_idx(ie, 2)) += H.col(1);
                    nodal_force.row(element_idx(ie, 0)) += (-H.col(0) - H.col(1));
                }

                else if ((MODEL == "CoRotated")) {
                    Eigen::Matrix2d H = CoRotatedFemStep(element_pos, D_m(ie), F, strain, init_area, LAM, MU);
                    nodal_force.row(element_idx(ie, 1)) += H.col(0);
                    nodal_force.row(element_idx(ie, 2)) += H.col(1);
                    nodal_force.row(element_idx(ie, 0)) += (-H.col(0) - H.col(1));
                }

                else if ((MODEL == "Neo")) {
                    Eigen::Matrix2d H = NeohookeanFemStep(element_pos, D_m(ie), F, strain, init_area, LAM, MU);
                    nodal_force.row(element_idx(ie, 1)) += H.col(0);
                    nodal_force.row(element_idx(ie, 2)) += H.col(1);
                    nodal_force.row(element_idx(ie, 0)) += (-H.col(0) - H.col(1));
                }
                else {
                    throw "Invalid Model";
                }

                if (TIME_INTE == "IMP") {
                    //ds, dF, dP, dH矩阵
                    set_zero_3dimtensor(dDsdx);
                    set_zero_3dimtensor(dF);
                    set_zero_3dimtensor(dP);
                    set_zero_3dimtensor(dH);

                    if (MODEL == "STVK") {
                        set_zero_3dimtensor(dE);
                        //D_s
                        calculate_D_s(dDsdx);
                        //dH
                        implicit_euler_stvk(D_m(ie), F, strain, dDsdx, dF, dE, dP, dH, init_area, LAM , MU);
                    }
                    else if (MODEL == "Neo") {
                        //D_s
                        calculate_D_s(dDsdx);
                        //dH
                        implicit_euler_Neo(D_m(ie), F, strain, dDsdx, dF, dP, dH, init_area, LAM, MU);
                    
                    }
                    else {
                        throw "Invalid Model";
                    }

                    //计算K矩阵
                    for (int n = 0; n < 3; n++) {
                        int node_idx = element_idx(ie, n);
                        for (int d = 0; d < 2; d++) {
                            int kidx = node_idx * 2 + d;
                            int didx = n * 2 + d;
                            int idx = element_idx(ie, 1) * 2;
                            K(idx, kidx) += dH(didx)(0,0);
                            //cout << K(idx, kidx) << endl;
                            idx = element_idx(ie, 1) * 2 + 1;
                            K(idx, kidx) += dH(didx)(1, 0);
                            //cout << K(idx, kidx) << endl;
                            idx = element_idx(ie, 2) * 2;
                            K(idx, kidx) += dH(didx)(0, 1);
                            //cout << K(idx, kidx) << endl;
                            idx = element_idx(ie, 2) * 2 + 1;
                            K(idx, kidx) += dH(didx)(1, 1);
                            //cout << K(idx, kidx) << endl;
                            idx = element_idx(ie, 0) * 2;
                            K(idx, kidx) += -dH(didx)(0, 0) - dH(didx)(0, 1);
                            //cout << K(idx, kidx) << endl;
                            idx = element_idx(ie, 0) * 2 + 1;
                            K(idx, kidx) += -dH(didx)( 1, 0) - dH(didx)(1, 1);
                            //cout << K(idx, kidx) << endl;
                        } 
                    }
                }
            }

            // 辛欧拉法更新位置
            if (TIME_INTE == "EXP") {
                for (int k = 0; k < NODE_NUM; k++) {
                    // Update position (Symplectic Euler)
                    node_vel.row(k) += (100 * (nodal_force.row(k) + gravity.row(0)) / POINT_MASS) * DT;
                    position.row(k) += DT * node_vel.row(k);
                    node_vel *= exp(-DT * 0.2);
                }
            }
            else if (TIME_INTE == "IMP") {
                // 加重力
                for (int node = 0; node < NODE_NUM; node++) {
                    nodal_force.row(node) += gravity.row(0);
                }
                nodal_force *= 100;
                
                // 计算X
                Eigen::MatrixXd X = solver(krow, NODE_NUM, POINT_MASS, DT, A, K, b, node_vel, nodal_force);
               
                for (int nod = 0; nod < NODE_NUM; nod++) {
                    node_vel(nod,0) = X(nod * 2 + 0);
                    node_vel(nod,1) = X(nod * 2 + 1);
                    position.row(nod) += node_vel.row(nod) * DT;
                    node_vel.row(nod) *= std::exp(-DT * 1.0);
                }
                
                K.setZero();
                A.setZero();
                b.setZero();
            }
            else
            {
                throw "Invalid Method";
            }

            // 边界条件
            for (int nodei = 0; nodei < NODE_NUM; nodei++) {
                if (position(nodei, 0) > 30) {
                    position(nodei, 0) = 30;
                    node_vel(nodei, 0) = -0.99 * node_vel(nodei, 0);
                }

                if (position(nodei, 1) > 15) {
                    position(nodei, 1) = 15;
                    node_vel(nodei, 1) = -0.99 * node_vel(nodei, 1);
                }

                if (position(nodei, 0) < 0) {
                    position(nodei, 0) = 0;
                    node_vel(nodei, 0) = -0.99 * node_vel(nodei, 0);
                }

                if (position(nodei, 1) < 0) {
                    position(nodei, 1) = 0;
                    node_vel(nodei, 1) = -0.99 * node_vel(nodei, 1);
                }

            }

        }
    }
    return 0;
}
