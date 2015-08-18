#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H
#include <vector>
#include "particle.h"

#define PI 3.141592f
#define INF 1E-12f
#define BOUNDARY 0.001f


class ParticleSystem {
public:
    ParticleSystem();
    ~ParticleSystem();
    void Step();
    float* Draw(float * res, float * colors);
    void createParticle();

protected:
    std::vector<Particle> Particles;
private:
    int column;
    int row;
    int num_particles;

    float GAS_STIFFNESS;
    float POLY6;
    float KERNEL;
    float TIMESTAMP;
    float REST_DENSITY;
    float SPIKY;
    float SELF_DENSITY;
    float GRAVITY_X;
    float GRAVITY_Y;
    float GRAD_POLY6;
    float MASS;
    float KERNEL2;
    float VISCOSITY;
    float pressure;
    int step_count;
    //float *res;


    float kernel(float x);
    float grad_kernel(float x);
};

#endif // PARTICLESYSTEM_H
