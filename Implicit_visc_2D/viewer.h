#ifndef VIEWER_H
#define VIEWER_H

#include <QGLWidget>
#include <QGLShaderProgram>
#include <iostream>
#include <QGLWidget>
#include <QGLShaderProgram>
#include <QMatrix4x4>
#include <QtGlobal>
#include <particlesystem.h>
#include <QTimer>

#if (QT_VERSION >= QT_VERSION_CHECK(5, 1, 0))
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#else
#include <QGLBuffer>
#endif
class Viewer : public QGLWidget {

    Q_OBJECT
public:
    Viewer(const QGLFormat& format, QWidget *parent = 0);
    virtual ~Viewer();
    QSize minimumSizeHint() const;
    QSize sizeHint() const;
signals:

public slots:
    void step();
    void createParticle();

protected:
    // Called when GL is first initialized
    virtual void initializeGL();
    // Called when our window needs to be redrawn
    virtual void paintGL();
    void resizeGL(int width, int height);
private:
    QMatrix4x4 getCameraMatrix();
    QMatrix4x4 mTransformMatrix;
    QTimer * mTimer;
    QTimer * mStepTimer;
    QTimer * createTimer;
    int MaxParticles;
    int ParticlesCount;
    GLuint billboard_vertex_buffer;
    GLuint particles_position_buffer;
    GLuint particles_color_buffer;
    QGLShaderProgram mProgram;
    GLuint vertexbuffer;
    GLint mvpMatrix;
    GLuint LoadShaders(const char * vertex_file_path,const char * fragment_file_path);
    QMatrix4x4 mPerspMatrix;
    ParticleSystem particleSystem;
#if (QT_VERSION >= QT_VERSION_CHECK(5, 1, 0))
    QOpenGLBuffer vbo;
    QOpenGLBuffer matbo;
    QOpenGLVertexArrayObject mVertexArrayObject;
#else
    QGLBuffer mVertexBufferObject;
#endif
};

#endif // VIEWER_H
