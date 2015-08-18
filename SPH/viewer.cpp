#include "viewer.h"
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <fstream>

Viewer::Viewer(const QGLFormat& format, QWidget *parent)
    : QGLWidget(format, parent),
      MaxParticles(2000),
      ParticlesCount(900)
#if (QT_VERSION >= QT_VERSION_CHECK(5, 1, 0))
    , vbo(QOpenGLBuffer::VertexBuffer)
    , matbo(QOpenGLBuffer::VertexBuffer)
    , mVertexArrayObject(this)
#else
    , vbo(QGLBuffer::VertexBuffer)
#endif
{
    mTimer = new QTimer(this);
    connect(mTimer, SIGNAL(timeout()), this, SLOT(update()));
    mTimer->start(1000/30);
    // Nothing to do here right now.

    mStepTimer = new QTimer(this);
    connect(mTimer, SIGNAL(timeout()), this, SLOT(step()));
    mTimer->start(1000/500);

    createTimer = new QTimer(this);
    connect(createTimer, SIGNAL(timeout()), this, SLOT(createParticle()));
    createTimer->start(1000/10);
}


void Viewer::createParticle() {
    if (ParticlesCount > 1200) {
        return;
    }
    particleSystem.createParticle();
    ParticlesCount++;
}

void Viewer::step() {
    particleSystem.Step();
}

Viewer::~Viewer() {
    // Nothing to do here right now.
}

void Viewer::initializeGL() {
    // Do some OpenGL setup
    QGLFormat glFormat = QGLWidget::format();
    if (!glFormat.sampleBuffers()) {
        std::cerr << "Could not enable sample buffers." << std::endl;
        return;
    }
    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);
    glClearColor(0.7, 0.7, 0.7, 0.0);


    int programID = LoadShaders("/Users/Aliya/SPH/shader.vert", "/Users/Aliya/SPH/shader.frag");
    glUseProgram(programID);
    static const GLfloat g_vertex_buffer_data[] = {
       -0.01f, -0.01f, 0.01f,
       0.01f, -0.01f, 0.01f,
       0.0f,  0.01f, 0.01f,
    };
    // This will identify our vertex buffer


    // Generate 1 buffer, put the resulting identifier in vertexbuffer
    glGenBuffers(1, &vertexbuffer);

    // The following commands will talk about our 'vertexbuffer' buffer
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);

    // Give our vertices to OpenGL.
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

    // The VBO containing the positions and sizes of the particles
    glGenBuffers(1, &particles_position_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
    // Initialize with empty (NULL) buffer : it will be updated later, each frame.
    glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);
    mvpMatrix = glGetUniformLocation(programID, "mvpMatrix");

    // The VBO containing the positions and sizes of the particles
    glGenBuffers(2, &particles_color_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
    // Initialize with empty (NULL) buffer : it will be updated later, each frame.
    glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);
}

void Viewer::paintGL() {
    //particleSystem.Step();
    glClearColor(0.7, 0.7, 0.7, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    //std::cerr << "no";
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(
       0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
       3,                  // size
       GL_FLOAT,           // type
       GL_FALSE,           // normalized?
       0,                  // stride
       (void*)0            // array buffer offset
    );


    // 2nd attribute buffer : positions of particles' centers
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
    glVertexAttribPointer(
     1, // attribute. No particular reason for 1, but must match the layout in the shader.
     4, // size : x + y + z + size => 4
     GL_FLOAT, // type
     GL_FALSE, // normalized?
     0, // stride
     (void*)0 // array buffer offset
    );
    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
    glVertexAttribPointer(
     2, // attribute. No particular reason for 1, but must match the layout in the shader.
     4, // size : x + y + z + size => 4
     GL_FLOAT, // type
     GL_FALSE, // normalized?
     0, // stride
     (void*)0 // array buffer offset
    );
    GLfloat g_particule_position_size_data[] = {
            0.0f, 0.0f, 0.0f, 0.5f,
            0.0f, 0.5f, 0.0f, 0.5f,
            0.0f, -0.5f, 0.0f, 0.5f,
            -0.5f, 0.0f, 0.0f, 0.5f,
            -0.5f, 0.5f, 0.0f, 0.5f,
            -0.5f, -0.5f, 0.0f, 0.5f,
            0.5f, 0.0f, 0.0f, 0.5f,
            0.5f, 0.5f, 0.0f, 0.5f,
            0.5f, -0.5f, 0.0f, 0.5f,
        };
    GLfloat g_particule_color_data[] = {
        1.0f, 0.0f, 0.0f
    };
    glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
    glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
    float particles[ParticlesCount*4];
    float colors[ParticlesCount*4];
    particleSystem.Draw(particles, colors);
    glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLfloat) * 4, particles);

    glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
    glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
    glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLfloat) * 4, colors);
    glVertexAttribDivisor(0, 0); // particles vertices : always reuse the same 4 vertices -> 0
    glVertexAttribDivisor(1, 1); // positions : one per quad (its center) -> 1
    glVertexAttribDivisor(2, 1);
    // Draw the triangle !
    //glDrawArrays(GL_TRIANGLES, 0, 3); // Starting from vertex 0; 3 vertices total -> 1 triangle

    QMatrix4x4 vMatrix;
    QMatrix4x4 cameraTransformation;
    QVector3D cameraPosition = cameraTransformation * QVector3D(0, 0, 1.0);
    QVector3D cameraUpDirection = cameraTransformation * QVector3D(0, 1, 0);
    QMatrix4x4  mvp = mPerspMatrix * vMatrix;
    vMatrix.lookAt(cameraPosition, QVector3D(0, 0, 0), cameraUpDirection);
    float vmmm[] = {
        1.7f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.7f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f,
    };
    //mProgram.setUniformValue(mvpMatrix, mPerspMatrix * vMatrix);
    glPointSize(3.0);
    //std::cerr << mvp.constData() << std::endl;
    glUniformMatrix4fv(mvpMatrix, 1, false, vmmm);

    glDrawArraysInstanced(GL_POINTS, 0, 1, ParticlesCount);
    //glDrawArrays(GL_TRIANGLE_STRIP, 0, 3);
    glDisableVertexAttribArray(0);
}


void Viewer::resizeGL(int width, int height) {
    if (height == 0) {
        height = 1;
    }

    mPerspMatrix.setToIdentity();
    mPerspMatrix.perspective(60.0, (float) width / (float) height, 0.001, 1000);

    glViewport(0, 0, width, height);
}

QSize Viewer::minimumSizeHint() const {
    return QSize(100, 100);
}

QSize Viewer::sizeHint() const {
    return QSize(600, 600);
}

GLuint Viewer::LoadShaders(const char * vertex_file_path,const char * fragment_file_path){

    // Create the shaders
    GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

    // Read the Vertex Shader code from the file
    std::string VertexShaderCode;
    std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
    if(VertexShaderStream.is_open())
    {
        std::string Line = "";
        while(getline(VertexShaderStream, Line))
            VertexShaderCode += "\n" + Line;
        VertexShaderStream.close();
    }

    // Read the Fragment Shader code from the file
    std::string FragmentShaderCode;
    std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
    if(FragmentShaderStream.is_open()){
        std::string Line = "";
        while(getline(FragmentShaderStream, Line))
            FragmentShaderCode += "\n" + Line;
        FragmentShaderStream.close();
    }

    GLint Result = GL_FALSE;
    int InfoLogLength;

    // Compile Vertex Shader
    printf("Compiling shader : %s\n", vertex_file_path);
    char const * VertexSourcePointer = VertexShaderCode.c_str();
    glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
    glCompileShader(VertexShaderID);

    // Check Vertex Shader
    glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    std::vector<char> VertexShaderErrorMessage(InfoLogLength);
    glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
    fprintf(stdout, "%s\n", &VertexShaderErrorMessage[0]);

    // Compile Fragment Shader
    printf("Compiling shader : %s\n", fragment_file_path);
    char const * FragmentSourcePointer = FragmentShaderCode.c_str();
    glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
    glCompileShader(FragmentShaderID);

    // Check Fragment Shader
    glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    std::vector<char> FragmentShaderErrorMessage(InfoLogLength);
    glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
    fprintf(stdout, "%s\n", &FragmentShaderErrorMessage[0]);

    // Link the program
    fprintf(stdout, "Linking program\n");
    GLuint ProgramID = glCreateProgram();
    glAttachShader(ProgramID, VertexShaderID);
    glAttachShader(ProgramID, FragmentShaderID);
    glLinkProgram(ProgramID);

    // Check the program
    glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
    glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    std::vector<char> ProgramErrorMessage( std::max(InfoLogLength, int(1)) );
    glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
    fprintf(stdout, "%s\n", &ProgramErrorMessage[0]);

    glDeleteShader(VertexShaderID);
    glDeleteShader(FragmentShaderID);

    return ProgramID;
}
