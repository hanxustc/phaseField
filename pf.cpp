#include "pf.h"
#include <cstdio>

PhaseField::PhaseField() {}
PhaseField1D::PhaseField1D() {}
PhaseField2D::PhaseField2D(ParameterPF* _parameterPF, ParameterInt* _parameterInt, ParameterMater* _parameterMater, InitialCondition* _initialCondition, MPIHandler* _mpiHandler) :
    parameterPF(_parameterPF), parameterInt(_parameterInt), parameterMater(_parameterMater), initialCondition(_initialCondition), mpiHandler(_mpiHandler) {}

void PhaseField2D::mpiSplit() {
    cols_per_process = parameterInt->Ny / mpiHandler->getSize();
    start_col = mpiHandler->getRank() * cols_per_process;
    end_col = start_col + cols_per_process;
    if (mpiHandler->getRank() == mpiHandler->getSize() - 1) {
        end_col = parameterInt->Ny;
    }
}

void PhaseField2D::allocate() {
    grad_theta.set_size(2 * parameterInt->Nx, parameterInt->Ny);
    norm2_gradtheta.set_size(parameterInt->Nx, parameterInt->Ny);
    norm_gradtheta.set_size(parameterInt->Nx, parameterInt->Ny);

    theta_smooth.set_size(parameterInt->Nx, parameterInt->Ny);
    theta_min.set_size(parameterInt->Nx, parameterInt->Ny);
    theta_max.set_size(parameterInt->Nx, parameterInt->Ny);
    dtheta.set_size(parameterInt->Nx, parameterInt->Ny);
    u1_new.set_size(parameterInt->Nx, parameterInt->Ny);
    alpha_new.set_size(parameterInt->Nx, parameterInt->Ny);
    Psi.set_size(parameterInt->Nx, parameterInt->Ny);

    w.set_size(parameterInt->Nx, parameterInt->Ny);
    dw.set_size(2 * parameterInt->Nx, parameterInt->Ny);
    B.set_size(parameterInt->Nx, parameterInt->Ny);
    b.set_size(parameterInt->Nx, parameterInt->Ny);
    db.set_size(2 * parameterInt->Nx, parameterInt->Ny);
    c.set_size(parameterInt->Nx, parameterInt->Ny);
    dc.set_size(parameterInt->Nx, parameterInt->Ny);
    dL_theta.set_size(parameterInt->Nx, parameterInt->Ny);
    dL_gradtheta.set_size(2 * parameterInt->Nx, parameterInt->Ny);
    df.set_size(parameterInt->Nx, parameterInt->Ny);
}

void PhaseField2D::smoothing() {
    mat theta_smooth_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            int m1 = i - 1; int n1 = jj - 1;
            int m2 = i + 1; int n2 = jj + 1;
            if (m1 < 0) {
                m1 = parameterInt->Nx - 1;
            }
            if (m2 > parameterInt->Nx - 1) {
                m2 = 0;
            }
            if (n1 < 0) {
                n1 = parameterInt->Ny - 1;
            }
            if (n2 > parameterInt->Ny - 1) {
                n2 = 0;
            }
            theta_smooth_local(i, j) = (parameterPF->theta(m1, jj) + parameterPF->theta(m2, jj) + parameterPF->theta(i, n1) + parameterPF->theta(i, n2) 
                                      + parameterPF->theta(m1, n1) + parameterPF->theta(m2, n1) + parameterPF->theta(m1, n2) + parameterPF->theta(m2, n2)
                                      + 8 * parameterPF->theta(i, jj)) / 16;
        }
    }
    MPI_Allgather(theta_smooth_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  theta_smooth.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
}

void PhaseField2D::smooth() {
    parameterPF->theta = theta_smooth;
}

void PhaseField2D::gradient() {
    mat grad_theta_local(2 * parameterInt->Nx, end_col - start_col, fill::zeros);
    mat norm2_gradtheta_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat norm_gradtheta_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            int m1 = i - 1; int n1 = jj - 1;
            int m2 = i + 1; int n2 = jj + 1;
            if (m1 < 0) {
                m1 = parameterInt->Nx - 1;
            }
            if (m2 > parameterInt->Nx - 1) {
                m2 = 0;
            }
            if (n1 < 0) {
                n1 = parameterInt->Ny - 1;
            }
            if (n2 > parameterInt->Ny - 1) {
                n2 = 0;
            }
            grad_theta_local(span(2 * i, 2 * i + 1), j) = { {(parameterPF->theta(m2, jj) - parameterPF->theta(m1, jj)) / parameterInt->dx / 2},
                                                            {(parameterPF->theta(i, n2) - parameterPF->theta(i, n1)) / parameterInt->dy / 2} };
            norm2_gradtheta_local(i, j) = pow(grad_theta_local(2 * i, j), 2) + pow(grad_theta_local(2 * i + 1, j), 2);
            norm_gradtheta_local(i, j) = pow(norm2_gradtheta_local(i, j), 0.5);
        }
    }
    MPI_Allgather(grad_theta_local.memptr(), 2 * parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  grad_theta.memptr(), 2 * parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(norm2_gradtheta_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  norm2_gradtheta.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(norm_gradtheta_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  norm_gradtheta.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
}

void PhaseField2D::inclination() {
    mat Psi_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            if (norm_gradtheta(i, jj) == 0) {
                Psi_local(i, j) = 0;
            } else if (grad_theta(2 * i, jj) == 0 && grad_theta(2 * i + 1, jj) > 0) {
                Psi_local(i, j) = pi / 2;
            } else if (grad_theta(2 * i, jj) == 0 && grad_theta(2 * i + 1, jj) < 0) {
                Psi_local(i, j) = -pi / 2;
            } else if (grad_theta(2 * i, jj) > 0) {
                Psi_local(i, j) = atan(grad_theta(2 * i + 1, jj) / grad_theta(2 * i, jj));
            } else if (grad_theta(2 * i, jj) < 0 && grad_theta(2 * i + 1, jj) < 0) {
                Psi_local(i, j) = -pi + atan(grad_theta(2 * i + 1, jj) / grad_theta(2 * i, jj));
            } else if (grad_theta(2 * i, jj) < 0 && grad_theta(2 * i + 1, jj) >= 0) {
                Psi_local(i, j) = pi + atan(grad_theta(2 * i + 1, jj) / grad_theta(2 * i, jj));
            }
        }
    }
    MPI_Allgather(Psi_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  Psi.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
}

void PhaseField2D::thetaPlusthetaMinusExtr2G() {
    mat theta_min_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat theta_max_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat dtheta_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            theta_min_local(i, j) = min({ initialCondition->theta_num(0), initialCondition->theta_num(1) });
            theta_max_local(i, j) = max({ initialCondition->theta_num(0), initialCondition->theta_num(1) });
            dtheta_local(i, j) = theta_max_local(i, j) - theta_min_local(i, j);
        }
    }
    MPI_Allgather(theta_min_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  theta_min.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(theta_max_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  theta_max.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(dtheta_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  dtheta.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
}

void PhaseField2D::thetaPlusthetaMinusExtr() {
    mat theta_min_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat theta_max_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat dtheta_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    
    int step = 6;
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            if (norm_gradtheta(i, jj) <= parameterMater->u0) {
                int min_index = abs(initialCondition->theta_num - parameterPF->theta(i, jj)).index_min();
                theta_min_local(i, j) = initialCondition->theta_num(min_index);
                theta_max_local(i, j) = initialCondition->theta_num(min_index);
            } else {
                mat v = { grad_theta(2 * i, jj) / norm_gradtheta(i, jj), grad_theta(2 * i + 1, jj) / norm_gradtheta(i, jj) };
                double dii1 = v(0), dii2 = -v(0), djj1 = v(1), djj2 = -v(1);
                double ii1 = 0, ii2 = 0, jj1 = 0, jj2 = 0;
                int iii1 = i, jjj1 = jj, iii2 = i, jjj2 = jj;
                
                while (norm_gradtheta(iii1, jjj1) > parameterMater->u0 && (abs(ii1) <= step && abs(jj1) <= step)) {
                    if (dii1 >= 0) {
                        iii1 = i + floor(ii1);
                    } else {
                        iii1 = i + ceil(ii1);
                    }
                    if (djj1 >= 0) {
                        jjj1 = jj + floor(jj1);
                    } else {
                        jjj1 = jj + ceil(jj1);
                    }
                    if (iii1 > parameterInt->Nx - 1) {
                        iii1 -= parameterInt->Nx;
                    } else if (iii1 < 0) {
                        iii1 += parameterInt->Nx;
                    }
                    if (jjj1 > parameterInt->Ny - 1) {
                        jjj1 -= parameterInt->Ny;
                    } else if (jjj1 < 0) {
                        jjj1 += parameterInt->Ny;
                    }
                    ii1 += dii1;
                    jj1 += djj1;
                }
                while (norm_gradtheta(iii2, jjj2) > parameterMater->u0 && (abs(ii2) <= step && abs(jj2) <= step)) {
                    if (dii2 >= 0) {
                        iii2 = i + floor(ii2);
                    } else {
                        iii2 = i + ceil(ii2);
                    }
                    if (djj2 >= 0) {
                        jjj2 = jj + floor(jj2);
                    } else {
                        jjj2 = jj + ceil(jj2);
                    }
                    if (iii2 > parameterInt->Nx - 1) {
                        iii2 -= parameterInt->Nx;
                    } else if (iii2 < 0) {
                        iii2 += parameterInt->Nx;
                    }
                    if (jjj2 > parameterInt->Ny - 1) {
                        jjj2 -= parameterInt->Ny;
                    } else if (jjj2 < 0) {
                        jjj2 += parameterInt->Ny;
                    }
                    ii2 += dii2;
                    jj2 += djj2;
                }
                int index1 = abs(initialCondition->theta_num - parameterPF->theta(iii1, jjj1)).index_min();
                int index2 = abs(initialCondition->theta_num - parameterPF->theta(iii2, jjj2)).index_min();
                theta_min_local(i, j) = min({ initialCondition->theta_num(index1), initialCondition->theta_num(index2) });
                theta_max_local(i, j) = max({ initialCondition->theta_num(index1), initialCondition->theta_num(index2) });
                dtheta_local(i, j) = theta_max_local(i, j) - theta_min_local(i, j);
            }
        }
    }
    MPI_Allgather(theta_min_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  theta_min.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(theta_max_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  theta_max.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(dtheta_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  dtheta.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
}

void PhaseField2D::newTerm2G() {
    mat alpha_new_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat u1_new_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            alpha_new_local(i, j) = parameterMater->alpha / pow(dtheta(i, jj), 2);
            u1_new_local(i, j) = parameterMater->u1 * dtheta(i, jj);
        }
    }
    MPI_Allgather(alpha_new_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  alpha_new.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(u1_new_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  u1_new.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
}

void PhaseField2D::newTerm() {
    mat alpha_new_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat u1_new_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            if (dtheta(i, jj) == 0) {
                alpha_new_local(i, j) = parameterMater->alpha;
                u1_new_local(i, j) = parameterMater->u1;
            } else {
                alpha_new_local(i, j) = parameterMater->alpha / pow(dtheta(i, jj), 2);
                u1_new_local(i, j) = parameterMater->u1 * dtheta(i, jj);
            }
        }
    }
    MPI_Allgather(alpha_new_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  alpha_new.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(u1_new_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  u1_new.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
}

void PhaseField2D::weightFunction() {
    mat w_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat dw_local(2 * parameterInt->Nx, end_col - start_col, fill::zeros);
    
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            if (norm_gradtheta(i, jj) <= parameterMater->u0) {
                w_local(i, j) = 0;
                dw_local(span(2 * i, 2 * i + 1), j) = { 0, 0 };
            } else if (norm_gradtheta(i, jj) >= u1_new(i ,jj)) {
                w_local(i, j) = 1;
                dw_local(span(2 * i, 2 * i + 1), j) = { 0, 0 };
            } else {
                w_local(i, j) = pow(norm_gradtheta(i, jj) - parameterMater->u0, 2) * pow(norm_gradtheta(i, jj) - 2 * u1_new(i, jj) + parameterMater->u0, 2) / pow(u1_new(i, jj) - parameterMater->u0, 4);
                dw_local(span(2 * i, 2 * i + 1), j) = 4 * (norm_gradtheta(i, jj) - parameterMater->u0) * (norm_gradtheta(i, jj) - 2 * u1_new(i, jj) + parameterMater->u0) * (norm_gradtheta(i, jj) - u1_new(i, jj)) / pow(u1_new(i, jj) - parameterMater->u0, 4) * grad_theta(span(2 * i, 2 * i + 1), jj) / norm_gradtheta(i, jj);
            }
        }
    }
    MPI_Allgather(w_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  w.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(dw_local.memptr(), 2 * parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  dw.memptr(), 2 * parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
}

void PhaseField2D::gbEnergy() {
    mat B_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat b_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat db_local(2 * parameterInt->Nx, end_col - start_col, fill::zeros);
    
    mat thetaw = (theta_min + theta_max) / 2 + parameterMater->K;
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            B_local(i, j) = abs(sin(parameterMater->n * dtheta(i, jj)));
            if (dtheta(i, jj) == 0) {
                b_local(i, j) = 0;
                db_local(span(2 * i, 2 * i + 1), j) = { 0, 0 };
            } else {
                b_local(i, j) = parameterMater->B0 * abs(sin(parameterMater->n * dtheta(i, jj))) / pow(dtheta(i, jj), 2) * (1 + parameterMater->epsilon * sin(parameterMater->m * (thetaw(i, jj) - Psi(i, jj))));
                if (grad_theta(2 * i, j) == 0) {
                    db_local(span(2 * i, 2 * i + 1), j) = { 0, 0 };
                } else {
                    db_local(span(2 * i, 2 * i + 1), j) = parameterMater->B0 * abs(sin(parameterMater->n * dtheta(i, jj))) / pow(dtheta(i, jj), 2) * parameterMater->m * parameterMater->epsilon * cos(parameterMater->m * (thetaw(i, jj) - Psi(i, jj))) / norm2_gradtheta(i, jj) * parameterMater->R * grad_theta(span(2 * i, 2 * i + 1), jj);
                }
            }
        }
    }
    MPI_Allgather(B_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  B.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(b_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  b.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(db_local.memptr(), 2 * parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  db.memptr(), 2 * parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
}

void PhaseField2D::doubleWell() {
    mat c_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat dc_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            if (dtheta(i, jj) == 0) {
                c_local(i, j) = 0;
                dc_local(i, j) = 0;
            } else {
                double u = (parameterPF->theta(i, jj) - theta_min(i, jj)) / dtheta(i, jj);
                c_local(i, j) = parameterMater->beta * pow(u, 2) * pow(1 - u, 2);
                dc_local(i, j) = parameterMater->beta * 2 * u * (1 - 3 * u + 2 * pow(u, 2)) / dtheta(i, jj);
            }
        }
    }
    MPI_Allgather(c_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  c.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(dc_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  dc.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
}

void PhaseField2D::dFunctional() {
    mat dL_theta_local(parameterInt->Nx, end_col - start_col, fill::zeros);
    mat dL_gradtheta_local(2 * parameterInt->Nx, end_col - start_col, fill::zeros);
    mat df_local(parameterInt->Nx, end_col - start_col, fill::zeros);

    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            dL_theta_local(i, j) = B(i, jj) * (1 - w(i, jj)) * dc(i, jj);
            dL_gradtheta_local(span(2 * i, 2 * i + 1), j) = b(i, jj) * (2 * w(i, jj) * grad_theta(span(2 * i, 2 * i + 1), jj) + dw(span(2 * i, 2 * i + 1), jj) * norm2_gradtheta(i, jj)) + w(i, jj) * db(span(2 * i, 2 * i + 1), jj) * norm2_gradtheta(i, jj)
                                                          + B(i, jj) * alpha_new(i, jj) * (2 * (1 - w(i, jj)) * grad_theta(span(2 * i, 2 * i + 1), jj) - dw(span(2 * i, 2 * i + 1), jj) * norm2_gradtheta(i, jj))
                                                          - B(i, jj) * dw(span(2 * i, 2 * i + 1), jj) * c(i, jj);
        }
    }
    MPI_Allgather(dL_theta_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  dL_theta.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(dL_gradtheta_local.memptr(), 2 * parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  dL_gradtheta.memptr(), 2 * parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
    
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < end_col - start_col; j++) {
            int jj = mpiHandler->getRank() * cols_per_process + j;
            
            int m1 = i - 1; int n1 = jj - 1;
            int m2 = i + 1; int n2 = jj + 1;
            if (m1 < 0) {
                m1 = parameterInt->Nx - 1;
            }
            if (m2 > parameterInt->Nx - 1) {
                m2 = 0;
            }
            if (n1 < 0) {
                n1 = parameterInt->Ny - 1;
            }
            if (n2 > parameterInt->Ny - 1) {
                n2 = 0;
            }
            df_local(i, j) = dL_theta(i, jj) - ((dL_gradtheta(2 * m2, jj) - dL_gradtheta(2 * m1, jj)) / (2 * parameterInt->dx) + (dL_gradtheta(2 * i + 1, n2) - dL_gradtheta(2 * i + 1, n1)) / (2 * parameterInt->dy));
        }
    }
    MPI_Allgather(df_local.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE,
                  df.memptr(), parameterInt->Nx * (end_col - start_col), MPI_DOUBLE, MPI_COMM_WORLD);
}

void PhaseField2D::totalEnergy() {
    F1 = 0; F2 = 0; F3 = 0;
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < parameterInt->Ny; j++) {
            F1 += b(i, j) * w(i, j) * norm2_gradtheta(i, j) * parameterInt->dx * parameterInt->dy;
            F2 += B(i, j) * alpha_new(i, j) * (1 - w(i, j)) * norm2_gradtheta(i, j) * parameterInt->dx * parameterInt->dy;
            F3 += B(i, j) * (1 - w(i, j)) * c(i, j) * parameterInt->dx * parameterInt->dy;
        }
    }
    F = F1 + F2 + F3;
}

void PhaseField2D::fdm() {
    allocate();
    mpiSplit();
    
    FILE* fid1 = fopen("theta.out", "w");;
    FILE* fid2 = fopen("theta12_num.out", "w");;
    FILE* fid3 = fopen("nodes.out", "w");
    FILE* fid4 = fopen("energy.out", "w");
    if (fid1 == nullptr) {
        cerr << "Failed to load file 1" << endl;
        return;
    }
    if (fid2 == nullptr) {
        cerr << "Failed to load file 2" << endl;
        return;
    }
    if (fid3 == nullptr) {
        cerr << "Failed to load file 3" << endl;
        return;
    }
    if (fid4 == nullptr) {
        cerr << "Failed to load file 4" << endl;
        return;
    }
    
    for (int i = 0; i < parameterInt->Nx; i++) {
        for (int j = 0; j < parameterInt->Ny; j++) {
            fprintf(fid1, "%10.6f\n", parameterPF->theta(i, j));
            fprintf(fid3, "%10.6f %10.6f\n", i + 0.5, j + 0.5);
        }
    }
    
    if (initialCondition->ngrain > 2) {
        for (int k = 1; k <= parameterInt->nstep; k++) {
            if (k % parameterInt->nsmooth == 0) {
                smoothing();
                smooth();
            }
            gradient();
            inclination();
            if (k % parameterInt->nthetapm == 0 || k == 1) {
                thetaPlusthetaMinusExtr();
            }
            newTerm();
            weightFunction();
            gbEnergy();
            doubleWell();
            dFunctional();
            parameterPF->theta += -parameterInt->dt * parameterMater->M * df;
            
            if (k % parameterInt->nprint == 0) {
                double F_temp = F;
                smoothing();
                totalEnergy();
                if (mpiHandler->getRank() == 0) {
                    for (int i = 0; i < parameterInt->Nx; i++) {
                        for (int j = 0; j < parameterInt->Ny; j++) {
                            fprintf(fid1, "%10.6f\n", theta_smooth(i, j));
                            fprintf(fid2, "%10.6f %10.6f\n", theta_min(i, j), theta_max(i, j));
                        }
                    }
                    fprintf(fid4, "%10.6f %10.6f %10.6f %10.6f\n", F1, F2, F3, F);
                    cout << k << endl;
                }
                
                /*When to stop*/
                if (abs(F - F_temp) < parameterInt->err) {
                    break;
                }
            }
        }
    } else if (initialCondition->ngrain == 2) {
        for (int k = 1; k <= parameterInt->nstep; k++) {
            if (k % parameterInt->nsmooth == 0) {
                smoothing();
                smooth();
            }
            gradient();
            inclination();
            if (k % parameterInt->nthetapm == 0 || k == 1) {
                thetaPlusthetaMinusExtr2G();
            }
            newTerm2G();
            weightFunction();
            gbEnergy();
            doubleWell();
            dFunctional();
            parameterPF->theta += -parameterInt->dt * parameterMater->M * df;
            
            if (k % parameterInt->nprint == 0) {
                double F_temp = F;
                smoothing();
                totalEnergy();
                if (mpiHandler->getRank() == 0) {
                    for (int i = 0; i < parameterInt->Nx; i++) {
                        for (int j = 0; j < parameterInt->Ny; j++) {
                            fprintf(fid1, "%10.6f\n", theta_smooth(i, j));
                        }
                    }
                    fprintf(fid2, "%10.6f %10.6f\n", theta_min(0, 0), theta_max(0, 0));
                    fprintf(fid4, "%10.6f %10.6f %10.6f %10.6f\n", F1, F2, F3, F);
                    cout << k << endl;
                }
                
                /*When to stop*/
                if (abs(F - F_temp) < parameterInt->err) {
                    break;
                }
            }
        }
    } else if (initialCondition->ngrain == 1) {
        if (mpiHandler->getRank() == 0) {
            fprintf(fid2, "%10.6f %10.6f\n", initialCondition->theta_num(0), initialCondition->theta_num(1));
            fprintf(fid4, "%10.6f %10.6f %10.6f %10.6f\n", 0, 0, 0, 0);
            cout << "No GB exists!" << endl;
        }
    }
    fclose(fid1);
    fclose(fid2);
    fclose(fid3);
    fclose(fid4);
}