#include "particlesystem.h"
#include <cmath>
#include <iostream>
using namespace std;
ParticleSystem::ParticleSystem()
{
    GAS_STIFFNESS = 1.0f  * 10;
    KERNEL = 0.04f;
    POLY6 = 315.0f/(64.0f * PI * pow(KERNEL, 9));
    TIMESTAMP = 0.003f;
    REST_DENSITY = 1000.0f;
    SPIKY = -45.0f/(PI * pow(KERNEL, 6));
    GRAVITY_X = 0.0f;
    GRAVITY_Y = 0.0f;
    GRAD_POLY6 = -945/(32 * PI * pow(KERNEL, 9));
    MASS = 0.018;
    KERNEL2 = KERNEL * KERNEL;
    SELF_DENSITY = MASS*POLY6*pow(KERNEL, 6);

    step_count = 0;
    VISCOSITY=6.5f/200.0f / 2.0f;

    row = 30;
    column = 30;
    num_particles = row*column;

    int layers = 15;
    num_particles = 4 * (std::pow(1.5f, layers) - 1);
    num_particles = layers * (2 * 3 + (layers - 1) * 3) / 2;


    for (int i = 0; i < layers; i++) {
        int num_per_layer = (i + 1) * 3;
        for (int j = 0; j < num_per_layer; j++) {
            float radius = KERNEL * 0.35f * float(i + 1.0f);
            float x = 0.0f + radius * cos(float(j) * 2.0f * 3.14f / float(num_per_layer));
            float y = 0.25f + radius * sin(float(j) * 2.0f * 3.14f / float(num_per_layer));
            cerr << i << " " << radius  << " " << num_per_layer << " "  << x << " " <<  y << endl;
            Particles.push_back(Particle(x, y, 0.5f, 0.0f, 0.2f, 0.7f));
        }

    }
}

void ParticleSystem::createParticle() {
    Particles.push_back(Particle(-0.3f + num_particles % 35 * 0.45f * KERNEL, 1.2f, 0.5f, 1.0f, 0.0f, 0.0f));
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

    if (step_count == 20) {
        GRAVITY_Y = -1.8f;
    }
    step_count++;
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
            Particles[i].r = 0.0f;
            Particles[i].g = 1.0f;
            Particles[i].b = 0.0f;
        }
    }

    // compute pressure force
    for (int i = 0; i < num_particles; i++) {
        F_i[i][0] = 0.0f;
        F_i[i][1] = 0.0f;
        for (int j = 0; j < num_particles; j++) {
            if (neighbors[i][j] < 0) continue;
            float r = sqrt(neighbors[i][j]);
            float temp_force;

            // VISCOSITY
            temp_force = MASS * VISCOSITY * 2 * MASS / rho[j]  * r / (r * r + 0.01 * KERNEL2) * 945/ (32 * PI * pow(KERNEL, 5)) * (-pow(1-r*r/KERNEL2, 2));
            F_i[i][0] += temp_force * (Particles[i].vx - Particles[j].vx);
            F_i[i][1] += temp_force * (Particles[i].vy - Particles[j].vy);
        }



    }

    // Predict advection velocity & Compute d_ii
    for (int i = 0; i < num_particles; i++) {
        Particles[i].vx += (F_i[i][0] + GRAVITY_X) * TIMESTAMP / MASS / 30.0f;
        Particles[i].vy += (F_i[i][1] + GRAVITY_Y * 2.0f) * TIMESTAMP / MASS /30.0f;
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
    for (int iter = 0; iter < 8; iter++) {
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

        if(Particles[i].y >= 0.5f)
        {
            Particles[i].vy = Particles[i].vy*wall_damping;
            Particles[i].y = 0.5f;
        }

        if(Particles[i].y < -0.5f)
        {
            Particles[i].vy =Particles[i].vy *wall_damping;
            Particles[i].y = -0.5f + BOUNDARY;
        }

        Particles[i].evx = (Particles[i].evx + Particles[i].vx) / 2;
        Particles[i].evy = (Particles[i].evy + Particles[i].vy) / 2;
    }
}
