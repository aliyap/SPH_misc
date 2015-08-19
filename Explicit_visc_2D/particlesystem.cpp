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
    TIMESTAMP = 0.0000009;
    REST_DENSITY = 1000.0f;
    GRAVITY_X = 0.0f;
    GRAVITY_Y = 0.0f;
    MASS = 0.018;
    KERNEL2 = KERNEL * KERNEL;
    SELF_DENSITY = MASS*POLY6*pow(KERNEL, 6);
    VISCOSITY= 6.5f/ 200.0f * 1700000.0f;//40.0f;

    row = 160;
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
            float y = 2.5f - i * 0.45 * KERNEL;
            Particles.push_back(Particle(x, y, 0.5f, r, g, b)); //, r * 200.0f - 100.0f, g * 200.0f - 100.0f));
        }

    }

    step_num = 0;

    mu_i = *(new vector<float>(num_particles));
    wij_x = *(new vector<float>(num_particles));
    alphaij_x = *(new vector<float>(num_particles));
    wij_y = *(new vector<float>(num_particles));
    alphaij_y = *(new vector<float>(num_particles));
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
        cerr << "hAHA" << endl;
    }
    step_num++;
    float VOLUME = 0.0001f;



    // find neighbors
    float F_i[num_particles][2];
    float rho[num_particles];
    float pressure_i[num_particles];
    float neighbors[num_particles][num_particles];

    float positions[num_particles][2];
    for (int i = 0; i < num_particles; i++) {
       int num_neighbors = 0;

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


    // compute rho
    for (int i = 0; i < num_particles; i++) {
        rho[i] = 0;

        for (int j = 0; j < num_particles; j++) {
            if (neighbors[i][j] < 0) continue;
            float r =  sqrt(neighbors[i][j]);
            rho[i] += MASS * 315.0f/ (64.0f *PI * pow(KERNEL, 3.0f)) * pow(1.0f-r*r/KERNEL2, 3.0f) ;
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

    }

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletListX;
    tripletListX.reserve(num_particles * 2);


    for (int i = 0; i < num_particles; i++) {
        Particles[i].vy += (GRAVITY_Y * 2.0f) * TIMESTAMP / MASS ;//30.0f;
        //cerr << i << " " << Particles[i].r << " " << old_vel(2 * i) <<  endl;
    }

    // EXPLICIT VISCOSITY
    // find s
    float S[num_particles][4];
    for (int i = 0; i < num_particles; i++) {
        S[i][0] = S[i][1] = S[i][2] = S[i][3] = 0.0f;
        for (int j = 0; j < num_particles; j++) {
            if ((neighbors[i][j]) < 0.0f) continue;
            float r =  sqrt(neighbors[i][j]);
            float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-r*r/KERNEL2, 2.0f));
            float delta_vx = (Particles[j].vx - Particles[i].vx);
            float delta_vy = Particles[j].vy - Particles[i].vy;
            float Wij_x = gradient * (Particles[i].x - Particles[j].x);
            float Wij_y = gradient * (Particles[i].y - Particles[j].y);
            float mu_v = VISCOSITY * MASS / rho[j];
            S[i][0] += mu_v * ( delta_vx * Wij_x + Wij_x * delta_vx);
            S[i][1] +=  mu_v * ( delta_vx * Wij_y + Wij_x * delta_vy);
            S[i][2] +=  mu_v * ( delta_vy * Wij_x + Wij_y * delta_vx);
            S[i][3] +=  mu_v * ( delta_vy * Wij_y + Wij_y * delta_vy);
            cerr << VISCOSITY << " " << MASS << " " << rho[j] << " " << mu_v << endl;
        }
    }

    // Compute velocity
    for (int i = 0; i < num_particles; i++) {
        float temp_x = 0.0f;
        float temp_y = 0.0f;
        for (int j = 0; j < num_particles; j++) {
            if ((neighbors[i][j]) < 0.0f) continue;
            float r = sqrt(neighbors[i][j]);
            float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-r*r/KERNEL2, 2.0f));
            float inner_sum[4] = {0.0f, 0.0f, 0.0f, 0.0f};
            float Wij_x = gradient * (Particles[i].x - Particles[j].x);
            float Wij_y = gradient * (Particles[i].y - Particles[j].y);
            for (int l = 0; l < 4; l++) {
                inner_sum[l] += S[i][l] / (rho[i] * rho[i]) + S[j][l] / (rho[j] * rho[j]);
            }
            temp_x += inner_sum[0] * Wij_x + inner_sum[1] * Wij_y;
            temp_y += inner_sum[2] * Wij_x + inner_sum[3] * Wij_y;
        }
        Particles[i].vx += MASS * TIMESTAMP * temp_x;
        Particles[i].vy += MASS * TIMESTAMP * temp_y;
    }



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

            for (int k = 0; k < num_particles; k++) {
                float C_i_grad_k[2]  = {0.0f, 0.0f};
                if (k == i) {
                    for (int j = 0; j < num_particles; j++) {
                        if ((neighbors[i][j]) < 0.0f) continue;
                        float r =  sqrt(neighbors[i][j]);
                        float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-r*r/KERNEL2, 2.0f));
                        C_i_grad_k[0] += gradient * (positions[i][0] - positions[j][0]) / REST_DENSITY;
                        C_i_grad_k[1] += gradient * (positions[i][1] - positions[j][1]) / REST_DENSITY;
                    }
                } else if (neighbors[i][k] > 0.0f){
                    float r =  sqrt(neighbors[i][k]);
                    float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-r*r/KERNEL2, 2.0f));
                    C_i_grad_k[0] = - gradient * (positions[i][0] - positions[k][0]) / REST_DENSITY;
                    C_i_grad_k[1] = - gradient * (positions[i][1] - positions[k][1]) / REST_DENSITY;
                }

                denom += pow( C_i_grad_k[0], 2) + pow( C_i_grad_k[1], 2);


            }
            lambda[i] = - C_i / (denom + 100.0f);
            //cerr << "denom " << i << "\t" << denom << endl;
        }

        float delta_p[num_particles][2];
        for (int i = 0; i < num_particles; i++) {
            delta_p[i][0] = 0.0f;
            delta_p[i][1] = 0.0f;
            for (int j = 0; j < num_particles; j++) {
                if (neighbors[i][j] < 0.0f) continue;
                float r =  sqrt(neighbors[i][j]);
                float gradient = 945.0f / (32.0f * PI * pow(KERNEL, 5.0f)) * (-pow(1.0f-r*r/KERNEL2, 2.0f));
                delta_p[i][0] += (lambda[i] + lambda[j]) * gradient * (positions[i][0] - positions[j][0]) / REST_DENSITY;
                delta_p[i][1] += (lambda[i] + lambda[j]) * gradient * (positions[i][1] - positions[j][1]) / REST_DENSITY;
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

        if(Particles[i].y >= 7.0f)
        {
            Particles[i].vy = Particles[i].vy*wall_damping;
            Particles[i].y = 7.0f;
        }

        if(Particles[i].y < -0.55f)
        {
            Particles[i].vy =Particles[i].vy *wall_damping;
            Particles[i].y = -0.55f + BOUNDARY;
        }

        Particles[i].evx = (Particles[i].evx + Particles[i].vx) / 2;
        Particles[i].evy = (Particles[i].evy + Particles[i].vy) / 2;
    }
}
