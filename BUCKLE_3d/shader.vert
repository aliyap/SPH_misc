#version 330

layout (location = 0) in vec3 vert;
layout (location = 1 ) in vec4 posn;
layout (location = 2) in vec4 color;

uniform mat4 mvpMatrix;

out vec4 fragColor;

void main()
{
    //gl_PointSize = 10.0;
    fragColor = color;
    gl_Position =  mvpMatrix * vec4(vert[0] + posn[0], vert[1] + posn[1], vert[2] + posn[2], 1.0);

    //gl_Position = mvpMatrix * vec4(vert, 1.0);
}
