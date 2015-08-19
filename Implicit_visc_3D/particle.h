#ifndef PARTICLE_H
#define PARTICLE_H

#include <QVector3D>
class Particle {

    QVector3D posn;
    float density;
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
    float evx;
    float evy;
    float evz;
    float pressure;
    Particle(float x, float y, float z, float r, float g, float b);
};


#endif // PARTICLE_H
