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
    //TIMESTAMP = 0.00001;
    TIMESTAMP = 0.00007;
    REST_DENSITY = 1000.0f;
    GRAVITY_X = 0.0f;
    GRAVITY_Y = 0.0f;
    MASS = 0.018;
    KERNEL2 = KERNEL * KERNEL;
    SELF_DENSITY = MASS*POLY6*pow(KERNEL, 6);
    //VISCOSITY= 3.25f * 17000.0f;//40.0f;
    VISCOSITY= 350000.0f;//40.0f;

    row = 600;
    column = 2;
    num_particles = row*column; //+ 40;

    for (int i = 0; i < row; i++) {
        float offset = 0.5 * i * KERNEL;
        for (int j = 0; j < column; j++) {
            int k = rand() % 2 ;
            float r = k == 0 ? 1.0f : 0.0;
            float g = k == 1 ? 1.0f : 0.0;
            float b = k == 2 ? 0.0f : 0.0;
            float x = 0.0f + j * 0.45 * KERNEL;
            float y = -0.5f + i * 0.45 * KERNEL;
            Particles.push_back(Particle(x, y, 0.5f, r, g, b)); //, r * 200.0f - 100.0f, g * 200.0f - 100.0f));
        }

    }

    coefs = new float*[num_particles * 2];
    for (int i = 0; i < 2 * num_particles; i++) {
        coefs[i] = new float[num_particles * 2];
    }

    /*
    for (int i = 0; i < row; i++) {
        float offset = 0.5 * i * KERNEL;
        for (int j = 0; j < column; j++) {
            int k = rand() % 2 ;
            float r = k == 0 ? 1.0f : 0.0;
            float g = k == 1 ? 1.0f : 0.0;
            float b = k == 2 ? 0.0f : 0.0;
            float x = -0.25f + offset + j * 0.45 * KERNEL;
            float y = 0.25f - i * 0.45 * KERNEL;
            Particles.push_back(Particle(x, y, 0.5f, r, g, b)); //, r * 200.0f - 100.0f, g * 200.0f - 100.0f));
        }

    }*/

    step_num = 0;

    mu_i = *(new vector<float>(num_particles));
    wij_x = *(new vector<float>(num_particles));
    alphaij_x = *(new vector<float>(num_particles));
    wij_y = *(new vector<float>(num_particles));
    alphaij_y = *(new vector<float>(num_particles));


    // initialize grid
    grid_rows = ceil((30.0 + 0.5) / KERNEL);
    grid_cols = ceil((6) / KERNEL);
    grid = new vector<int>[grid_rows * grid_cols];
}

int ParticleSystem::offset(int row, int column) {
    return row * grid_cols + column;
}

int ParticleSystem::posn_to_grid(float x, float y, float z) {
    int r = floor((y + 0.5f) / KERNEL);
    int c = floor((x + 3.0f) / KERNEL);
    return offset(r, c);
}

int ParticleSystem::posn_to_grid(float x, float y, float z, int &r, int &c) {
    r = floor((y + 0.5f) / KERNEL);
    c = floor((x + 3.0f) / KERNEL);
    return 0;
}

void ParticleSystem::createParticle() {
    Particles.push_back(Particle(0.0f + num_particles % 5 * 0.45f * KERNEL, 1.2f, 0.5f, 1.0f, 0.0f, 0.0f));
    num_particles++;
}

ParticleSystem::~ParticleSystem() {
    //delete res;
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
    float neighbors[num_particles][num_particles];

    float positions[num_particles][2];
    for (int i = 0; i < num_particles; i++) {
        for (int j = 0; j < num_particles; j++) {
            neighbors[i][j] = -1.0f;
            if (i == j) {
                continue;
            }
            float d = pow(Particles[i].x - Particles[j].x, 2) + pow(Particles[i].y - Particles[j].y, 2);
            if (d < KERNEL2) {;
                neighbors[i][j] = d;
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
        int r, c;
        posn_to_grid(Particles[i].x, Particles[i].y, Particles[i].z, r, c);

        for (int row = r - 1; row <= r + 1; row++) {
            for (int col = c - 1; col <= c + 1; col++) {
                if (row < 0 || row >= grid_rows || col < 0 && col >= grid_cols) break;
                int index = offset(row, col);

                for (int p = 0; p < grid[index].size(); p++) {
                    int j = grid[index][p];
                    if (neighbors[i][j] < 0) continue;
                    float r =  sqrt(neighbors[i][j]);
                    rho[i] += MASS * 315.0f/ (64.0f *PI * pow(KERNEL, 3.0f)) * pow(1.0f-r*r/KERNEL2, 3.0f) ;
                }
            }
        }
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
        //rho[i] = 1000.0f;
    }


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

    for (int i = 0; i < num_particles * 2; i++) {
        for (int j = 0; j < num_particles * 2; j++) {
            coefs[i][j] = 0.0f;
        }
    }

    float mass_t = - MASS * TIMESTAMP;

    // Compute velocity
    for (int i = 0; i < num_particles; i++) {
        float sum_wij_x = 0.0f;
        float sum_wij_y = 0.0f;

        int r, c;
        posn_to_grid(Particles[i].x, Particles[i].y, Particles[i].z, r, c);

        for (int row = r - 1; row <= r + 1; row++) {
            for (int col = c - 1; col <= c + 1; col++) {
                if (row < 0 || row >= grid_rows || col < 0 && col >= grid_cols) break;
                int index = offset(row, col);

                for (int p = 0; p < grid[index].size(); p++) {
                    int j = grid[index][p];
                    if ((neighbors[i][j]) < 0.0f) continue;
                    float rdist =  sqrt(neighbors[i][j]);
                    float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-rdist*rdist/KERNEL2, 2.0f));
                    float Wij_x = gradient * (Particles[i].x - Particles[j].x);
                    float Wij_y = gradient * (Particles[i].y - Particles[j].y);
                    sum_wij_x += Wij_x;
                    sum_wij_y += Wij_y;
                }
            }
        }

        for (int row = r - 1; row <= r + 1; row++) {
            for (int col = c - 1; col <= c + 1; col++) {
                if (row < 0 || row >= grid_rows || col < 0 && col >= grid_cols) break;
                int index = offset(row, col);

                for (int p = 0; p < grid[index].size(); p++) {
                    int j = grid[index][p];
                    if ((neighbors[i][j]) < 0.0f) continue;
                    float rdist =  sqrt(neighbors[i][j]);
                    float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-rdist*rdist/KERNEL2, 2.0f));
                    float Wij_x = gradient * (Particles[i].x - Particles[j].x);
                    float Wij_y = gradient * (Particles[i].y - Particles[j].y);
                    float mu_v = VISCOSITY * MASS / rho[j];
                    //Si[0] += mu_v/ (rho[i] * rho[i]) * ( 2.0f * Particles[j].vx * Wij_x - 2.0f * Particles[i].vx * Wij_x) * sum_wij_x;
                    //Si[1] +=  mu_v/ (rho[i] * rho[i]) * ( Particles[j].vx * Wij_y - Particles[i].vx * Wij_y + Particles[j].vy * Wij_x - Particles[i].vy * Wij_x) * sum_wij_y;

                    coefs[i * 2][j * 2] += mass_t * mu_v / (rho[i] * rho[i]) * (2.0f * Wij_x * sum_wij_x + Wij_y * sum_wij_y);
                    coefs[i * 2][i * 2] -= mass_t * mu_v/ (rho[i] * rho[i]) * (2.0f * Wij_x * sum_wij_x + Wij_y * sum_wij_y);
                    coefs[i * 2][j * 2 + 1] += mass_t * mu_v/ (rho[i] * rho[i]) * Wij_x * sum_wij_y;
                    coefs[i * 2][i * 2 + 1] -= mass_t * mu_v/ (rho[i] * rho[i]) * Wij_x * sum_wij_y;

                    //Si[2] +=  mu_v/ (rho[i] * rho[i]) * ( Particles[j].vy * Wij_x - Particles[i].vy * Wij_x + Particles[j].vx * Wij_y - Particles[i].vx * Wij_y) * sum_wij_x;
                    //Si[3] +=  mu_v/ (rho[i] * rho[i]) * (2.0f *  Particles[j].vy * Wij_y - 2.0f * Particles[i].vy * Wij_y) * sum_wij_y;

                    coefs[i * 2 + 1][j * 2 + 1] += mass_t * mu_v/ (rho[i] * rho[i]) * (Wij_x * sum_wij_x + 2.0f * Wij_y * sum_wij_y);
                    coefs[i * 2 + 1][i * 2 + 1] -= mass_t * mu_v/ (rho[i] * rho[i]) * (Wij_x * sum_wij_x + 2.0f * Wij_y * sum_wij_y);
                    coefs[i * 2 + 1][j * 2] += mass_t * mu_v/ (rho[i] * rho[i]) * Wij_y * sum_wij_x;
                    coefs[i * 2 + 1][i * 2] -= mass_t * mu_v/ (rho[i] * rho[i]) * Wij_y * sum_wij_x;

                    int rk, ck;
                    posn_to_grid(Particles[j].x, Particles[j].y, Particles[j].z, rk, ck);

                    for (int rowk = rk - 1; rowk <= rk + 1; rowk++) {
                        for (int colk = ck - 1; colk <= ck + 1; colk++) {
                            if (rowk < 0 || rowk >= grid_rows || colk < 0 && colk >= grid_cols) break;
                            int indexk = offset(rowk, colk);

                            for (int pk = 0; pk < grid[indexk].size(); pk++) {
                                int k = grid[indexk][pk];
                                if ((neighbors[j][k]) < 0.0f) continue;
                                float rjk =  sqrt(neighbors[j][k]);
                                float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-rjk*rjk/KERNEL2, 2.0f));
                                float Wjk_x = gradient * (Particles[j].x - Particles[k].x);
                                float Wjk_y = gradient * (Particles[j].y - Particles[k].y);
                                float mu_v = VISCOSITY * MASS / rho[k];
                                //Sj[0] += mu_v / (rho[j] * rho[j]) * ( 2.0f * Particles[k].vx * Wjk_x - 2.0f * Particles[j].vx * Wjk_x)* Wij_x;
                                //Sj[1] +=  mu_v / (rho[j] * rho[j]) * ( Particles[k].vx * Wjk_y - Particles[j].vx * Wjk_y + Particles[k].vy * Wjk_x - Particles[j].vy * Wjk_x)* Wij_y;

                                coefs[i * 2][k * 2] += mass_t * mu_v / (rho[j] * rho[j]) * (2.0f * Wjk_x * Wij_x + Wjk_y * Wij_y);
                                coefs[i * 2][j * 2] -= mass_t * mu_v / (rho[j] * rho[j]) * (2.0f * Wjk_x * Wij_x + Wjk_y * Wij_y);
                                coefs[i * 2][k * 2 + 1] += mass_t * mu_v / (rho[j] * rho[j]) * Wjk_x * Wij_y;
                                coefs[i * 2][j * 2 + 1] -= mass_t * mu_v / (rho[j] * rho[j]) * Wjk_x * Wij_y;

                                //Sj[2] +=  mu_v / (rho[j] * rho[j]) * ( Particles[k].vy * Wjk_x - Particles[j].vy * Wjk_x + Particles[k].vx * Wjk_y - Particles[j].vx * Wjk_y)* Wij_x;
                                //Sj[3] +=  mu_v / (rho[j] * rho[j]) * (2.0f *  Particles[k].vy * Wjk_y - 2.0f * Particles[j].vy * Wjk_y)* Wij_y;

                                coefs[i * 2 + 1][k * 2 + 1] += mass_t * mu_v / (rho[j] * rho[j]) * (Wjk_x * Wij_x + 2.0f * Wjk_y * Wij_y);
                                coefs[i * 2 + 1][j * 2 + 1] -= mass_t * mu_v / (rho[j] * rho[j]) * (Wjk_x * Wij_x + 2.0f * Wjk_y * Wij_y);
                                coefs[i * 2 + 1][k * 2] += mass_t * mu_v / (rho[j] * rho[j]) * Wjk_y * Wij_x;
                                coefs[i * 2 + 1][j * 2] -= mass_t * mu_v / (rho[j] * rho[j]) * Wjk_y * Wij_x;
                            }
                        }
                    }
                }
            }
        }

    }


    for (int i = 0; i < num_particles; i++) {
        coefs[2 * i][2 * i] += 1.0f;
        coefs[2 * i + 1][2 * i + 1] += 1.0f;
    }


    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletListX;

    for (int i = 0; i < 2 * num_particles; i++) {
        for (int j = 0; j < 2 * num_particles; j++) {
            if (coefs[i][j] != 0.0f) {
                tripletListX.push_back(T(i, j, coefs[i][j]));
            }
        }
    }

    SparseMatrix<double, RowMajor> AX = SparseMatrix<double, RowMajor>(2 * num_particles, 2 * num_particles);
    // fill A
    VectorXd old_vel(2 * num_particles);
    for (int i = 0; i < num_particles; i++) {
      old_vel(2 * i) = Particles[i].vx;
      old_vel(2 * i + 1) = Particles[i].vy;
    }


    VectorXd new_vel(num_particles * 2);
    // fill b
    // solve Ax = b
    SparseLU<SparseMatrix<double, RowMajor> > solver;
    AX.setFromTriplets(tripletListX.begin(), tripletListX.end());
    solver.compute(AX);


    new_vel = solver.solve(old_vel);

    // Predict advection velocity & Compute d_ii
    for (int i = 0; i < num_particles; i++) {
        Particles[i].vx = new_vel(2 * i);
        Particles[i].vy = new_vel(2 * i + 1);
    }

    /*
    for (int i = 0; i < num_particles; i++) {
        for (int j = 0; j < num_particles; j++) {
            if (fabs(coefs[i][j] - coefs[j][i]) > 0.00001f) {
                cerr << i << "\t" << j << "\t" << coefs[i][j] << "\t" << coefs[j][i] << endl;
            }
        }
    }
    */

    /*for (int i = 0; i < num_particles; i++) {
        float temp_v = 0.0f;
        for (int j = 0; j < num_particles; j++) {
            temp_v += coefs[2 * i][2 * j] * Particles[j].vx + coefs[2 * i][2 * j + 1] * Particles[j].vy;
        }
        cerr << temp_v << " " << velocities[i][0] << endl;
    }*/

/*
    for (int i = 0; i < num_particles; i++) {
        float temp_x = 0.0f;
        float temp_y = 0.0f;
        for (int l = 0; l < num_particles; l++) {
            temp_x += coefs[i * 2][2 * l] * velocities[l][0] + coefs[i * 2][2 * l + 1] * velocities[l][1];
            temp_y += coefs[i * 2 + 1][2 * l] * velocities[l][0] + coefs[i * 2 + 1][2 * l + 1] * velocities[l][1];
        }

        Particles[i].vx = temp_x; //Particles[i].vx +
        Particles[i].vy = temp_y; //Particles[i].vy +
    }
*/


    for (int i = 0; i < num_particles; i++) {
        positions[i][0] = Particles[i].x + Particles[i].vx * TIMESTAMP;
        positions[i][1] = Particles[i].y + Particles[i].vy * TIMESTAMP;
    }

    // Mark neighbours
    for (int i = 0; i < num_particles; i++) {
        for (int j = 0; j < num_particles; j++) {
            neighbors[i][j] = -1.0f;
            if (i == j) {
                continue;
            }
            float d = pow(positions[i][0] - positions[j][0], 2) + pow(positions[i][1] - positions[j][1], 2);
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

            float C_i_grad_k[2]  = {0.0f, 0.0f};

            int r, c;
            posn_to_grid(Particles[i].x, Particles[i].y, Particles[i].z, r, c);

            for (int row = r - 1; row <= r + 1; row++) {
                for (int col = c - 1; col <= c + 1; col++) {
                    if (row < 0 || row >= grid_rows || col < 0 && col >= grid_cols) break;
                    int index = offset(row, col);

                    for (int p = 0; p < grid[index].size(); p++) {
                        int j = grid[index][p];
                        if ((neighbors[i][j]) < 0.0f) continue;
                        float r =  sqrt(neighbors[i][j]);
                        float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-r*r/KERNEL2, 2.0f));
                        C_i_grad_k[0] += gradient * (positions[i][0] - positions[j][0]) / REST_DENSITY;
                        C_i_grad_k[1] += gradient * (positions[i][1] - positions[j][1]) / REST_DENSITY;
                    }
                }
            }
            denom += pow( C_i_grad_k[0], 2) + pow( C_i_grad_k[1], 2);

            posn_to_grid(Particles[i].x, Particles[i].y, Particles[i].z, r, c);

            for (int row = r - 1; row <= r + 1; row++) {
                for (int col = c - 1; col <= c + 1; col++) {
                    if (row < 0 || row >= grid_rows || col < 0 && col >= grid_cols) break;
                    int index = offset(row, col);

                    for (int p = 0; p < grid[index].size(); p++) {
                        int k = grid[index][p];
                        float C_i_grad_k[2]  = {0.0f, 0.0f};
                        if (neighbors[i][k] > 0.0f){
                            float r =  sqrt(neighbors[i][k]);
                            float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-r*r/KERNEL2, 2.0f));
                            C_i_grad_k[0] = - gradient * (positions[i][0] - positions[k][0]) / REST_DENSITY;
                            C_i_grad_k[1] = - gradient * (positions[i][1] - positions[k][1]) / REST_DENSITY;
                        }
                        denom += pow( C_i_grad_k[0], 2) + pow( C_i_grad_k[1], 2);
                    }
                }
            }
            lambda[i] = - C_i / (denom + 100.0f);
        }

        float delta_p[num_particles][2];
        for (int i = 0; i < num_particles; i++) {
            delta_p[i][0] = 0.0f;
            delta_p[i][1] = 0.0f;
            int r, c;
            posn_to_grid(Particles[i].x, Particles[i].y, Particles[i].z, r, c);

            for (int row = r - 1; row <= r + 1; row++) {
                for (int col = c - 1; col <= c + 1; col++) {
                    if (row < 0 || row >= grid_rows || col < 0 && col >= grid_cols) break;
                    int index = offset(row, col);

                    for (int p = 0; p < grid[index].size(); p++) {
                        int j = grid[index][p];
                        if (neighbors[i][j] < 0.0f) continue;
                        float r =  sqrt(neighbors[i][j]);
                        float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-r*r/KERNEL2, 2.0f));
                        delta_p[i][0] += (lambda[i] + lambda[j]) * gradient * (positions[i][0] - positions[j][0]) / REST_DENSITY;
                        delta_p[i][1] += (lambda[i] + lambda[j]) * gradient * (positions[i][1] - positions[j][1]) / REST_DENSITY;
                    }
                }
            }

            positions[i][0] += delta_p[i][0];
            positions[i][1] += delta_p[i][1];
        }

    }


    // change velocity and position
    for (int i = 0; i < num_particles; i++) {
        Particles[i].vx = (positions[i][0] - Particles[i].x) / TIMESTAMP;
        Particles[i].vy = (positions[i][1] - Particles[i].y) / TIMESTAMP;
        Particles[i].x = positions[i][0];
        Particles[i].y = positions[i][1];


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

        Particles[i].evx = (Particles[i].evx + Particles[i].vx) / 2;
        Particles[i].evy = (Particles[i].evy + Particles[i].vy) / 2;
    }

    for (int i = 0; i < grid_rows * grid_cols; i++) {
        grid[i].clear();
    }
}
