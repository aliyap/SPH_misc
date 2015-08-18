#include "particlesystem.h"
#include <cmath>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/SparseLU"
#include <math.h>
using namespace std;
using namespace Eigen;
ParticleSystem::ParticleSystem()
{
    KERNEL = 0.04f;
    TIMESTAMP = 0.00007;
    REST_DENSITY = 1000.0f;
    GRAVITY_X = 0.0f;
    GRAVITY_Y = 0.0f;
    MASS = 0.018;
    KERNEL2 = KERNEL * KERNEL;
    VISCOSITY= 400000.0f;

    row = 200;
    column = 2;
    // the third dimension (z)
    int thickness = 2;
    num_particles = thickness * row*column;

    for (int zc = -1; zc < 1; zc++) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                int k = rand() % 2 ;
                float x = 0.0f + j * 0.45 * KERNEL;
                float y = -0.5f + i * 0.55 * KERNEL;
                float z = 0.0 + zc * 0.35 * KERNEL;
                Particles.push_back(Particle(x, y, z, 0, 0, 0));
            }
        }
    }

    // coefficients of the sparse matrix
    coefs = new float*[num_particles * 3];
    // positions of the particles
    positions = new float*[num_particles];
    neighbors = new float*[num_particles];
    for (int i = 0; i < 3 * num_particles; i++) {
        coefs[i] = new float[num_particles * 3];
    }
    for (int i = 0; i < num_particles; i++) {
        positions[i] = new float[3];
        neighbors[i] = new float[num_particles];
    }
    step_num = 0;

    mu_i = *(new vector<float>(num_particles));
    wij_x = *(new vector<float>(num_particles));
    alphaij_x = *(new vector<float>(num_particles));
    wij_y = *(new vector<float>(num_particles));
    alphaij_y = *(new vector<float>(num_particles));

    // initialize grid
    grid_rows = ceil((30.0 + 0.5) / KERNEL);
    grid_cols = ceil((6) / KERNEL);
    grid_z = ceil(4.0f / KERNEL);
    grid = new vector<int>[grid_rows * grid_cols * grid_z];
}

int ParticleSystem::offset(int row, int column, int z) {
    return z * grid_cols * grid_rows + row * grid_cols + column; //
}

int ParticleSystem::posn_to_grid(float x, float y, float z) {
    int r = floor((y + 0.5f) / KERNEL);
    int c = floor((x + 3.0f) / KERNEL);
    int zc = floor((z + 2.0f) / KERNEL);
    return offset(r, c, zc);
}

int ParticleSystem::posn_to_grid(float x, float y, float z, int &r, int &c, int &zc) {
    r = floor((y + 0.5f) / KERNEL);
    c = floor((x + 3.0f) / KERNEL);
    zc = floor((z + 2.0f) / KERNEL);
    return 0;
}

void ParticleSystem::createParticle() {
    Particles.push_back(Particle(0.0f + num_particles % 5 * 0.45f * KERNEL, 1.2f, 0.5f, 1.0f, 0.0f, 0.0f));
    num_particles++;
}

ParticleSystem::~ParticleSystem() {
}

float * ParticleSystem::Draw(float * res, float * colors) {
    //float res[num_particles*4];
    for (int i = 0; i < num_particles; i++) {
        res[i*4] = Particles[i].x;
        res[i*4+1] = Particles[i].y;
        res[i*4+2] = Particles[i].z;
        res[i*4+3] = 1.0f;
        colors[i*4] = Particles[i].r;
        colors[i*4+1] = Particles[i].g;
        colors[i*4+2] = Particles[i].b;
        colors[i*4+3] = 1.0f;
    }
    return res;
}

void ParticleSystem::Step() {

    if (step_num == 2) {
        GRAVITY_Y = -1.8f * 100.0f;
    }
    step_num++;


    // find neighbors
    float rho[num_particles];

    for (int i = 0; i < num_particles; i++) {
        for (int j = 0; j < num_particles; j++) {
            neighbors[i][j] = -1.0f;
            if (i == j) {
                continue;
            }
            float d2 = pow(Particles[i].x - Particles[j].x, 2) + pow(Particles[i].y - Particles[j].y, 2) + pow(Particles[i].z - Particles[j].z, 2);
            if (d2 < KERNEL2) {;
                neighbors[i][j] = d2;
            }
        }
    }
    // need to assign to grid cells
    for (int i = 0; i < num_particles; i++) {
        int grid_index = posn_to_grid(Particles[i].x, Particles[i].y, Particles[i].z);
        grid[grid_index].push_back(i);
    }

    // compute rho
    for (int i = 0; i < num_particles; i++) {
        rho[i] = 0;
        int r, c, zc;
        posn_to_grid(Particles[i].x, Particles[i].y, Particles[i].z, r, c, zc);

        for (int row = r - 1; row <= r + 1; row++) {
            for (int col = c - 1; col <= c + 1; col++) {
                for (int zlocal = zc - 1; zlocal <= zc + 1; zlocal++) {
                    if (row < 0 || row >= grid_rows || col < 0 || col >= grid_cols || zlocal < 0 || zlocal >= grid_z) continue;
                    int index = offset(row, col, zlocal);

                    for (int p = 0; p < grid[index].size(); p++) {
                        int j = grid[index][p];
                        if (neighbors[i][j] < 0) continue;
                        float r =  sqrt(neighbors[i][j]);
                        rho[i] += MASS * 315.0f/ (64.0f *PI * pow(KERNEL, 3.0f)) * pow(1.0f-r*r/KERNEL2, 3.0f) ;
                    }
                }
            }
        }

        // colour of the particle based on the deviation from rest density
        if (rho[i] > REST_DENSITY + 50.0f ) {
            Particles[i].r = 1.0f;
            Particles[i].g = 0.0f;
            Particles[i].b = 0.0f;
        } else if (rho[i] < REST_DENSITY - 50.0f ) {
            Particles[i].r = 0.0f;
            Particles[i].g = 0.0f;
            Particles[i].b = 1.0f;
        } else {
            Particles[i].r = 0.05f;
            Particles[i].g = 0.5f;
            Particles[i].b = 0.05f;
        }
    }

    // I. External Forces
    for (int i = 0; i < num_particles; i++) {
        Particles[i].vy += (GRAVITY_Y * 2.0f) * TIMESTAMP / MASS;
    }

    // EXPLICIT VISCOSITY
    // find s
    /*float velocities[num_particles][2];
    for (int i = 0; i < num_particles; i++) {
        velocities[i][0] = Particles[i].vx;
        velocities[i][1] = Particles[i].vy;
    }*/

    for (int i = 0; i < num_particles * 3; i++) {
        for (int j = 0; j < num_particles * 3; j++) {
            coefs[i][j] = 0.0f;
        }
    }

    float mass_t = - MASS * TIMESTAMP;

    // Compute velocity

    for (int i = 0; i < num_particles; i++) {
        float sum_wij_x = 0.0f;
        float sum_wij_y = 0.0f;
        float sum_wij_z = 0.0f;

        int r, c, zc;
        posn_to_grid(Particles[i].x, Particles[i].y, Particles[i].z, r, c, zc);

        for (int row = r - 1; row <= r + 1; row++) {
            for (int col = c - 1; col <= c + 1; col++) {
                for (int zlocal = zc - 1; zlocal <= zc + 1; zlocal++) {
                    if (row < 0 || row >= grid_rows || col < 0 || col >= grid_cols || zlocal < 0 || zlocal >= grid_z) continue;
                    int index = offset(row, col, zlocal);

                    for (int p = 0; p < grid[index].size(); p++) {
                        int j = grid[index][p];
                        if ((neighbors[i][j]) < 0.0f) continue;
                        float rdist =  sqrt(neighbors[i][j]);
                        float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-rdist*rdist/KERNEL2, 2.0f));
                        float Wij_x = gradient * (Particles[i].x - Particles[j].x);
                        float Wij_y = gradient * (Particles[i].y - Particles[j].y);
                        float Wij_z = gradient * (Particles[i].z - Particles[j].z);
                        sum_wij_x += Wij_x;
                        sum_wij_y += Wij_y;
                        sum_wij_z += Wij_z;
                    }
                }
            }
        }

        for (int row = r - 1; row <= r + 1; row++) {
            for (int col = c - 1; col <= c + 1; col++) {
                for (int zlocal = zc - 1; zlocal <= zc + 1; zlocal++) {
                    if (row < 0 || row >= grid_rows || col < 0 || col >= grid_cols || zlocal < 0 || zlocal >= grid_z) continue;
                    int index = offset(row, col, zlocal);

                    for (int p = 0; p < grid[index].size(); p++) {
                        int j = grid[index][p];
                        if ((neighbors[i][j]) < 0.0f) continue;
                        float rdist =  sqrt(neighbors[i][j]);
                        float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-rdist*rdist/KERNEL2, 2.0f));
                        float Wij_x = gradient * (Particles[i].x - Particles[j].x);
                        float Wij_y = gradient * (Particles[i].y - Particles[j].y);
                        float Wij_z = gradient * (Particles[i].z - Particles[j].z);
                        float mu_v = VISCOSITY * MASS / rho[j];
                        // Coefficients for u_i
                        // C_uj_ui
                        coefs[i * 3][j * 3] += mass_t * mu_v / (rho[i] * rho[i]) * (2.0f * Wij_x * sum_wij_x + Wij_y * sum_wij_y + Wij_z * sum_wij_z);
                        // C_ui_ui
                        coefs[i * 3][i * 3] -= mass_t * mu_v/ (rho[i] * rho[i]) * (2.0f * Wij_x * sum_wij_x + Wij_y * sum_wij_y + Wij_z * sum_wij_z);
                        // C_vj_ui
                        coefs[i * 3][j * 3 + 1] += mass_t * mu_v/ (rho[i] * rho[i]) * Wij_x * sum_wij_y;
                        // C_vi_ui
                        coefs[i * 3][i * 3 + 1] -= mass_t * mu_v/ (rho[i] * rho[i]) * Wij_x * sum_wij_y;
                        coefs[i * 3][j * 3 + 2] += mass_t * mu_v/ (rho[i] * rho[i]) * Wij_x * sum_wij_z;
                        coefs[i * 3][i * 3 + 2] -= mass_t * mu_v/ (rho[i] * rho[i]) * Wij_x * sum_wij_z;

                        // Coefficients for v_i
                        // C_vj_vi
                        coefs[i * 3 + 1][j * 3 + 1] += mass_t * mu_v/ (rho[i] * rho[i]) * (Wij_x * sum_wij_x + 2.0f * Wij_y * sum_wij_y + Wij_z * sum_wij_z);
                        // C_vi_vi
                        coefs[i * 3 + 1][i * 3 + 1] -= mass_t * mu_v/ (rho[i] * rho[i]) * (Wij_x * sum_wij_x + 2.0f * Wij_y * sum_wij_y + Wij_z * sum_wij_z);
                        // C_uj_vi
                        coefs[i * 3 + 1][j * 3] += mass_t * mu_v/ (rho[i] * rho[i]) * Wij_y * sum_wij_x;
                        // C_ui_vi
                        coefs[i * 3 + 1][i * 3] -= mass_t * mu_v/ (rho[i] * rho[i]) * Wij_y * sum_wij_x;
                        coefs[i * 3 + 1][j * 3 + 2] += mass_t * mu_v/ (rho[i] * rho[i]) * Wij_y * sum_wij_z;
                        coefs[i * 3 + 1][i * 3 + 2] -= mass_t * mu_v/ (rho[i] * rho[i]) * Wij_y * sum_wij_z;

                        // Coefficients for w_i
                        coefs[i * 3 + 2][i * 3] -= mass_t * mu_v/ (rho[i] * rho[i]) * Wij_z * sum_wij_x;
                        coefs[i * 3 + 2][i * 3 + 1] -= mass_t * mu_v/ (rho[i] * rho[i]) * Wij_z * sum_wij_y;
                        coefs[i * 3 + 2][i * 3 + 2] -= mass_t * mu_v/ (rho[i] * rho[i]) * (Wij_x * sum_wij_x + Wij_y * sum_wij_y + 2.0f * Wij_z * sum_wij_z);
                        coefs[i * 3 + 2][j * 3] += mass_t * mu_v/ (rho[i] * rho[i]) * Wij_z * sum_wij_x;
                        coefs[i * 3 + 2][j * 3 + 1] += mass_t * mu_v/ (rho[i] * rho[i]) * Wij_z * sum_wij_y;
                        coefs[i * 3 + 2][j * 3 + 2] += mass_t * mu_v/ (rho[i] * rho[i]) * (Wij_x * sum_wij_x + Wij_y * sum_wij_y + 2.0f * Wij_z * sum_wij_z);

                        int rk, ck, zck;
                        posn_to_grid(Particles[j].x, Particles[j].y, Particles[j].z, rk, ck, zck);

                        // for all neighbours k of j
                        for (int rowk = rk - 1; rowk <= rk + 1; rowk++) {
                            for (int colk = ck - 1; colk <= ck + 1; colk++) {
                                for (int zlocalk = zck - 1; zlocalk <= zck + 1; zlocalk++) {
                                if (rowk < 0 || rowk >= grid_rows || colk < 0 || colk >= grid_cols || zlocalk < 0 || zlocalk >= grid_z) continue;
                                int indexk = offset(rowk, colk, zlocalk);

                                    for (int pk = 0; pk < grid[indexk].size(); pk++) {
                                        int k = grid[indexk][pk];
                                        if ((neighbors[j][k]) < 0.0f) continue;
                                        float rjk =  sqrt(neighbors[j][k]);
                                        float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-rjk*rjk/KERNEL2, 2.0f));
                                        float Wjk_x = gradient * (Particles[j].x - Particles[k].x);
                                        float Wjk_y = gradient * (Particles[j].y - Particles[k].y);
                                        float Wjk_z = gradient * (Particles[j].z - Particles[k].z);
                                        float mu_v = VISCOSITY * MASS / rho[k];

                                        // Coefficients for u_i
                                        // C_uk_ui
                                        coefs[i * 3][k * 3] += mass_t * mu_v / (rho[j] * rho[j]) * (2.0f * Wjk_x * Wij_x + Wjk_y * Wij_y + Wjk_z * Wij_z);
                                        // C_uj_ui
                                        coefs[i * 3][j * 3] -= mass_t * mu_v / (rho[j] * rho[j]) * (2.0f * Wjk_x * Wij_x + Wjk_y * Wij_y + Wjk_z * Wij_z);
                                        // C_vk_ui
                                        coefs[i * 3][k * 3 + 1] += mass_t * mu_v / (rho[j] * rho[j]) * Wjk_x * Wij_y;
                                        // C_vj_ui
                                        coefs[i * 3][j * 3 + 1] -= mass_t * mu_v / (rho[j] * rho[j]) * Wjk_x * Wij_y;
                                        coefs[i * 3][k * 3 + 2] += mass_t * mu_v / (rho[j] * rho[j]) * Wjk_x * Wij_z;
                                        coefs[i * 3][j * 3 + 2] -= mass_t * mu_v / (rho[j] * rho[j]) * Wjk_x * Wij_z;

                                        // Coefficients for v_i
                                        coefs[i * 3 + 1][j * 3] -= mass_t * mu_v / (rho[j] * rho[j]) * Wjk_y * Wij_x;
                                        coefs[i * 3 + 1][j * 3 + 1] -= mass_t * mu_v / (rho[j] * rho[j]) * (Wjk_x * Wij_x + 2.0f * Wjk_y * Wij_y + Wjk_z * Wij_z);
                                        coefs[i * 3 + 1][j * 3 + 2] -= mass_t * mu_v / (rho[j] * rho[j]) * Wjk_y * Wij_z;
                                        coefs[i * 3 + 1][k * 3] += mass_t * mu_v / (rho[j] * rho[j]) * Wjk_y * Wij_x;
                                        coefs[i * 3 + 1][k * 3 + 1] += mass_t * mu_v / (rho[j] * rho[j]) * (Wjk_x * Wij_x + 2.0f * Wjk_y * Wij_y + Wjk_z * Wij_z);
                                        coefs[i * 3 + 1][k * 3 + 2] += mass_t * mu_v / (rho[j] * rho[j]) * Wjk_y * Wij_z;

                                        // Coefficients for w_i
                                        coefs[i * 3 + 2][j * 3] -= mass_t * mu_v / (rho[j] * rho[j]) * Wjk_z * Wij_x;
                                        coefs[i * 3 + 2][j * 3 + 1] -= mass_t * mu_v / (rho[j] * rho[j]) * Wjk_z * Wij_y;
                                        coefs[i * 3 + 2][j * 3 + 2] -= mass_t * mu_v / (rho[j] * rho[j]) * (Wjk_x * Wij_x + Wjk_y * Wij_y + 2.0f * Wjk_z * Wij_z);
                                        coefs[i * 3 + 2][k * 3] += mass_t * mu_v / (rho[j] * rho[j]) * Wjk_z * Wij_x;
                                        coefs[i * 3 + 2][k * 3 + 1] += mass_t * mu_v / (rho[j] * rho[j]) * Wjk_z * Wij_y;
                                        coefs[i * 3 + 2][k * 3 + 2] += mass_t * mu_v / (rho[j] * rho[j]) * (Wjk_x * Wij_x + Wjk_y * Wij_y + 2.0f * Wjk_z * Wij_z);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }


    // diagonal entries: TODO: iterating is very inefficient, do it some place else
    for (int i = 0; i < num_particles; i++) {
        coefs[3 * i][3 * i] += 1.0f;
        coefs[3 * i + 1][3 * i + 1] += 1.0f;
        coefs[3 * i + 2][3 * i + 2] += 1.0f;
    }

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletListX;

    for (int i = 0; i < 3 * num_particles; i++) {
        for (int j = 0; j < 3 * num_particles; j++) {
            if (coefs[i][j] != 0.0f) {
                tripletListX.push_back(T(i, j, coefs[i][j]));
            }
        }
    }

    SparseMatrix<double, RowMajor> AX = SparseMatrix<double, RowMajor>(3 * num_particles, 3 * num_particles);
    // fill A
    VectorXd old_vel(3 * num_particles);
    for (int i = 0; i < num_particles; i++) {
      old_vel(3 * i) = Particles[i].vx;
      old_vel(3 * i + 1) = Particles[i].vy;
      old_vel(3 * i + 2) = Particles[i].vz;
    }


    VectorXd new_vel(num_particles * 3);
    // fill b
    // solve Ax = b
    SparseLU<SparseMatrix<double, RowMajor> > solver;
    AX.setFromTriplets(tripletListX.begin(), tripletListX.end());
    solver.compute(AX);

    new_vel = solver.solve(old_vel);

    for (int i = 0; i < num_particles; i++) {
        Particles[i].vx = new_vel(3 * i);
        Particles[i].vy = new_vel(3 * i + 1);
        Particles[i].vz = new_vel(3 * i + 2);
    }

    for (int i = 0; i < num_particles; i++) {
        positions[i][0] = Particles[i].x + Particles[i].vx * TIMESTAMP;
        positions[i][1] = Particles[i].y + Particles[i].vy * TIMESTAMP;
        positions[i][2] = Particles[i].z + Particles[i].vz * TIMESTAMP;
    }

    // Mark neighbours
    for (int i = 0; i < num_particles; i++) {
        for (int j = 0; j < num_particles; j++) {
            neighbors[i][j] = -1.0f;
            if (i == j) {
                continue;
            }
            float d = pow(positions[i][0] - positions[j][0], 2) + pow(positions[i][1] - positions[j][1], 2) + pow(positions[i][2] - positions[j][2], 2);
            if (d < KERNEL2) {;
                neighbors[i][j] = d;
            }
        }
    }

    float lambda[num_particles];
    for (int iter = 0; iter < 7; iter++) {
        for (int i = 0; i < num_particles; i++) {
            float C_i = rho[i] / REST_DENSITY - 1.0f;
            float denom = 0.0f;

            float C_i_grad_k[3]  = {0.0f, 0.0f, 0.0f};

            int r, c, zc;
            posn_to_grid(Particles[i].x, Particles[i].y, Particles[i].z, r, c, zc);
            for (int row = r - 1; row <= r + 1; row++) {
                for (int col = c - 1; col <= c + 1; col++) {
                    for (int zlocal = zc - 1; zlocal <= zc + 1; zlocal++) {
                        if (row < 0 || row >= grid_rows || col < 0 || col >= grid_cols || zlocal < 0 || zlocal >= grid_z) continue;
                        int index = offset(row, col, zlocal);

                        for (int p = 0; p < grid[index].size(); p++) {
                            int j = grid[index][p];
                            if ((neighbors[i][j]) < 0.0f) continue;
                            float r =  sqrt(neighbors[i][j]);
                            float gradient = 945.0f / (64.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-r*r/KERNEL2, 2.0f));
                            C_i_grad_k[0] += gradient * (positions[i][0] - positions[j][0]) / REST_DENSITY;
                            C_i_grad_k[1] += gradient * (positions[i][1] - positions[j][1]) / REST_DENSITY;
                            C_i_grad_k[2] += gradient * (positions[i][2] - positions[j][2]) / REST_DENSITY;
                        }
                    }
                }
            }
            denom += pow( C_i_grad_k[0], 2) + pow( C_i_grad_k[1], 2) + pow( C_i_grad_k[2], 2);

            posn_to_grid(Particles[i].x, Particles[i].y, Particles[i].z, r, c, zc);

            for (int row = r - 1; row <= r + 1; row++) {
                for (int col = c - 1; col <= c + 1; col++) {
                    for (int zlocal = zc - 1; zlocal <= zc + 1; zlocal++) {
                        if (row < 0 || row >= grid_rows || col < 0 || col >= grid_cols || zlocal < 0 || zlocal >= grid_z) continue;
                        int index = offset(row, col, zlocal);

                        for (int p = 0; p < grid[index].size(); p++) {
                            int k = grid[index][p];
                            float C_i_grad_k[2]  = {0.0f, 0.0f};
                            if (neighbors[i][k] > 0.0f){
                                float r =  sqrt(neighbors[i][k]);
                                float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-r*r/KERNEL2, 2.0f));
                                C_i_grad_k[0] = - gradient * (positions[i][0] - positions[k][0]) / REST_DENSITY;
                                C_i_grad_k[1] = - gradient * (positions[i][1] - positions[k][1]) / REST_DENSITY;
                                C_i_grad_k[2] = - gradient * (positions[i][2] - positions[k][2]) / REST_DENSITY;
                            }
                            denom += pow( C_i_grad_k[0], 2) + pow( C_i_grad_k[1], 2) + pow( C_i_grad_k[2], 2);
                        }
                    }
                }
            }
            lambda[i] = - C_i / (denom + 100.0f);
        }

        float delta_p[num_particles][3];
        for (int i = 0; i < num_particles; i++) {
            delta_p[i][0] = 0.0f;
            delta_p[i][1] = 0.0f;
            delta_p[i][2] = 0.0f;
            int r, c, zc;
            posn_to_grid(Particles[i].x, Particles[i].y, Particles[i].z, r, c, zc);

            for (int row = r - 1; row <= r + 1; row++) {
                for (int col = c - 1; col <= c + 1; col++) {
                    for (int zlocal = zc - 1; zlocal <= zc + 1; zlocal++) {
                        if (row < 0 || row >= grid_rows || col < 0 || col >= grid_cols || zlocal < 0 || zlocal >= grid_z) continue;
                        int index = offset(row, col, zlocal);

                        for (int p = 0; p < grid[index].size(); p++) {
                            int j = grid[index][p];
                            if (neighbors[i][j] < 0.0f) continue;
                            float r =  sqrt(neighbors[i][j]);
                            float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-r*r/KERNEL2, 2.0f));
                            delta_p[i][0] += (lambda[i] + lambda[j]) * gradient * (positions[i][0] - positions[j][0]) / REST_DENSITY;
                            delta_p[i][1] += (lambda[i] + lambda[j]) * gradient * (positions[i][1] - positions[j][1]) / REST_DENSITY;
                            delta_p[i][2] += (lambda[i] + lambda[j]) * gradient * (positions[i][2] - positions[j][2]) / REST_DENSITY;
                        }
                    }
                }
            }

            positions[i][0] += delta_p[i][0];
            positions[i][1] += delta_p[i][1];
            positions[i][2] += delta_p[i][2];
        }

    }


    // change velocity and position
    for (int i = 0; i < num_particles; i++) {
        Particles[i].vx = (positions[i][0] - Particles[i].x) / TIMESTAMP;
        Particles[i].vy = (positions[i][1] - Particles[i].y) / TIMESTAMP;
        //cerr << positions[i][2] << " " << Particles[i].z << endl;
        Particles[i].vz = (positions[i][2] - Particles[i].z) / TIMESTAMP;
        Particles[i].x = positions[i][0];
        Particles[i].y = positions[i][1];
        Particles[i].z = positions[i][2];


        if (Particles[i].y < -7.0f)  Particles[i].y = -7.0f;

        float wall_damping=-0.5f;
        if(Particles[i].x >= 0.3f)
        {
            Particles[i].vx = Particles[i].vx * wall_damping;
            Particles[i].x = 0.3f - BOUNDARY;
        }

        if(Particles[i].x < -0.3f)
        {
            Particles[i].vx = Particles[i].x*wall_damping;
            Particles[i].x= -0.3f + BOUNDARY;
        }

        if(Particles[i].y >= 30.0f)
        {
            Particles[i].vy = Particles[i].vy*wall_damping;
            Particles[i].y = 30.0f;
        }

        if(Particles[i].y < -0.5f)
        {
            Particles[i].vy =Particles[i].vy *wall_damping;
            Particles[i].y = -0.5f + BOUNDARY;
        }

        if(Particles[i].z >= 2.0f)
        {
            Particles[i].vz = Particles[i].vz*wall_damping;
            Particles[i].z = 2.0f;
        }

        if(Particles[i].z < -2.0f)
        {
            Particles[i].vz =Particles[i].vz *wall_damping;
            Particles[i].z = -2.0f + BOUNDARY;
        }
    }

    for (int i = 0; i < grid_rows * grid_cols * grid_z; i++) {
        grid[i].clear();
    }
}
