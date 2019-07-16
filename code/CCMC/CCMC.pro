TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    jelliumbasis.cpp \
    cartesianstate.cpp

HEADERS += \
    abstractbasis.h \
    jelliumbasis.h \
    cartesianstate.h
