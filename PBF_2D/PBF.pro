#-------------------------------------------------
#
# Project created by QtCreator 2015-05-25T17:24:44
#
#-------------------------------------------------

QT       += opengl core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = PBF
TEMPLATE = app


SOURCES += main.cpp\
        appwindow.cpp \
    particle.cpp \
    particlesystem.cpp \
    viewer.cpp

HEADERS  += appwindow.h \
    particlesystem.h \
    particle.h \
    viewer.h

FORMS    += appwindow.ui

DISTFILES += \
    shader.frag \
    shader.vert
