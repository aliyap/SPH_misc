#include "particlesystem.h"
#include <cmath>
#include <iostream>
using namespace std;
ParticleSystem::ParticleSystem()
{
    GAS_STIFFNESS = 1.0f  * 10;
    KERNEL = 0.04f;
    POLY6 = 315.0f/(64.0f * PI * pow(KERNEL, 3));
    TIMESTAMP = 0.003f;
    REST_DENSITY = 1000.0f;
    SPIKY = -45.0f/(PI * pow(KERNEL, 6));
    GRAVITY_X = 0.0f;
    GRAVITY_Y = -1.8f;
    GRAD_POLY6 = 945/ (32 * PI * pow(KERNEL, 5));
    MASS = 0.018;
    KERNEL2 = KERNEL * KERNEL;
    SELF_DENSITY = MASS*POLY6*pow(KERNEL, 6);

    VISCOSITY=6.5f/ 200.0f;

    row = 30;
    column = 30;
    num_particles = row*column;

    // place the initial particles
    for (int i = 0; i < column; i++) {
        for (int j = 0; j < row; j++) {
            float x = -0.25f + KERNEL*0.45f * i;
            float y = 0.0f + KERNEL * 0.45f * j;
            float r,g, b;
            if (i < 15 && j < 15){
                r = 0.0f; g = 0.0f; b = 1.0f;
            } else if (i < 15 && j >= 15) {
                r = 1.0f; g = 0.0f; b = 1.0f;
            } else if (i >= 15 && j < 15) {
                r = 1.0f; g = 1.0f; b = 0.0f;
            } else {
                r = 0.0f; g = 1.0f; b = 1.0f;
            }
            Particles.push_back(Particle(x, y, 0.5f, r, g, b));
        }
    }

}

void ParticleSystem::createParticle() {
    Particles.push_back(Particle(-0.3f + num_particles % 35 * 0.45f * KERNEL, 1.2f, 0.5f, 1.0f, 0.0f, 0.0f));
    num_particles++;
}

ParticleSystem::~ParticleSystem() {
}

float * ParticleSystem::Draw(float * res, float * colors) {
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

float ParticleSystem::kernel(float x, float h) {
    float q = abs(x)/h;
    if (q >= 1) return 0;
    return 315.0f / (64.0f * 3.14f * pow(h, 3.0f)) * pow(1.0f-pow(q, 2.0f), 3.0f);
}

float ParticleSystem::grad_kernel(float x, float h) {
    float q = abs(x)/h;
    if (q>= 1) return 0;
    return 945.0f / (64.0f * 3.14f * pow(h, 5.0f)) * x * (-pow(1-q * q, 2.0f));
}

void ParticleSystem::Step() {
    // find neighbors
    float F_i[num_particles][2];
    float rho[num_particles];
    float pressure_i[num_particles];
    float neighbors[num_particles][num_particles];
    for (int i = 0; i < num_particles; i++) {
        F_i[i][0] = 0;
        F_i[i][1] = 0;
        rho[i] = 0;
        pressure_i[i] = 0;

        for (int j = 0; j < num_particles; j++) {
            neighbors[i][j] = -1.0f;
        }
    }


    for (int i = 0; i < num_particles; i++) {
        for (int j = 0; j < num_particles; j++) {
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
            rho[i] += MASS * POLY6 * pow(1-r*r/KERNEL2, 3) ;
        }

        // colour based on the density
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

        // advection force
        // VISCOSITY - based on Ihmsen et al, SPH in Computer Graphics, algorithm 1

        /*
        float temp_x = 0;
        float temp_y = 0;
        for (int j = 0; j < num_particles; j++) {
            if (neighbors[i][j] < 0) continue;
            float r = sqrt(neighbors[i][j]);
            float grad_kernel_x = GRAD_POLY6 * (Particles[i].vx - Particles[j].vx) * (- pow(1 - r * r / KERNEL2, 2));
            float grad_kernel_y = GRAD_POLY6 * (Particles[i].vy - Particles[j].vy) * (- pow(1 - r * r / KERNEL2, 2));
            temp_x += 2 * MASS / rho[j] * (Particles[i].vx - Particles[j].vx) * r * grad_kernel_x / (r * r + 0.01 * KERNEL2);
            temp_x += 2 * MASS / rho[j] * (Particles[i].vx - Particles[j].vx) * r * grad_kernel_x / (r * r + 0.01 * KERNEL2);
        }
        F_i[i][0] += MASS * VISCOSITY * temp_x;
        F_i[i][1] += MASS * VISCOSITY * temp_y;

        F_i[i][0] += temp_force * (Particles[i].vx - Particles[j].vx) ;
        F_i[i][1] += temp_force * (Particles[i].vy - Particles[j].vy) ;
        Particles[i].vx += (F_i[i][0] + GRAVITY_X) * TIMESTAMP / MASS ;
        Particles[i].x += Particles[i].vx * TIMESTAMP;
        Particles[i].vy += (F_i[i][1] + (-9.8)/90.0f ) * TIMESTAMP / MASS ;
        */


        pressure_i[i] = GAS_STIFFNESS * (pow(rho[i] / REST_DENSITY, 7) -1);
    }

    // compute pressure force
    for (int i = 0; i < num_particles; i++) {
        F_i[i][0] = 0.0f;
        F_i[i][1] = 0.0f;
        for (int j = 0; j < num_particles; j++) {
            if (neighbors[i][j] < 0) continue;
            //float v = MASS / rho[i] * 2;
            float r = sqrt(neighbors[i][j]);
            //float pres_kernel = SPIKY * pow(KERNEL - r, 2);
            float temp_force; //= v  * (pressure_i[i]/ pow(rho[i], 2) + pressure_i[j]/ pow(rho[j], 2)) * pres_kernel;
            temp_force = MASS * MASS * (pressure_i[i]/ pow(rho[i], 2) + pressure_i[j]/ pow(rho[j], 2)) * 45.0f / (PI * pow(KERNEL, 5)) * (- pow(1-r / KERNEL, 2)/ (r/ KERNEL));
            F_i[i][0] -= (Particles[i].x - Particles[j].x ) * temp_force ;// r;
            F_i[i][1] -= (Particles[i].y - Particles[j].y) * temp_force ;// r;

            // VISCOSITY
            float dist_x = (Particles[i].x - Particles[j].x);
            float dist_y = (Particles[i].y - Particles[j].y);
            float temp_force_x = MASS * VISCOSITY * MASS / rho[j] * 2  * dist_x * dist_x / (dist_x * dist_x + 0.01 * KERNEL2) * GRAD_POLY6 * (-pow(1-r*r/KERNEL2, 2));
            float temp_force_y = MASS * VISCOSITY * MASS / rho[j] * 2  * dist_y * dist_y / (dist_y * dist_y + 0.01 * KERNEL2) * GRAD_POLY6 * (-pow(1-r*r/KERNEL2, 2));

            F_i[i][0] += temp_force_x * (Particles[i].vx - Particles[j].vx) ;
            F_i[i][1] += temp_force_y * (Particles[i].vy - Particles[j].vy) ;

        }

    }

    // change velocity and position
    for (int i = 0; i < num_particles; i++) {
        Particles[i].vx += (F_i[i][0] + GRAVITY_X) * TIMESTAMP / MASS ;
        Particles[i].x += Particles[i].vx * TIMESTAMP;
        Particles[i].vy += (F_i[i][1] + (-9.8)/90.0f ) * TIMESTAMP / MASS ;
        //cerr << "Velocity " << Particles[i].vy << endl;
        Particles[i].y += Particles[i].vy * TIMESTAMP;
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

    }
}
