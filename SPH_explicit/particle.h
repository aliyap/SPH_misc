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
    Particle(float x, float y, float z, float r, float g, float b, float vx = 0.0f, float vy = 0.0f, float vz = 0.0f);
};


#endif // PARTICLE_H
