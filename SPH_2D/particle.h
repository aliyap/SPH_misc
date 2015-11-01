#ifndef PARTICLE_H
#define PARTICLE_H

#include <QVector3D>
class Particle {
public:
    float x;
    float y;
    float z;
    float r;
    float g;
    float b;
    float vx;
    float vy;
    float vz;
    Particle(float x, float y, float z, float r, float g, float b);
};


#endif // PARTICLE_H
