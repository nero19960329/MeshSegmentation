#include "Utils.h"

#include <cmath>

double distanceBetweenVectors(unsigned char a[4], unsigned char b[4]) {
    double tmp = 0.0;
    for (int i = 0; i < 4; ++i) {
        tmp += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sqrt(tmp);
}

bool equals(unsigned char a[4], unsigned char b[4], double k) {
    return ((unsigned char)(b[0] * k) == a[0]) && 
        ((unsigned char)(b[1] * k) == a[1]) && 
        ((unsigned char)(b[2] * k) == a[2]);
}

unsigned char* HSVtoRGB(double h, double s, double v) {
    h *= 360.0;

    int tmp = floor(h / 60);
    double f = h / 60 - tmp;
    double p = v * (1 - s);
    double q = v * (1 - f * s);
    double t = v * (1 - (1 - f) * s);

    double* tmpArray = new double[3];
    unsigned char* res = new unsigned char[4];

    if (tmp == 0) {
        tmpArray[0] = v;
        tmpArray[1] = t;
        tmpArray[2] = p;
    } else if (tmp == 1) {
        tmpArray[0] = q;
        tmpArray[1] = v;
        tmpArray[2] = p;
    } else if (tmp == 2) {
        tmpArray[0] = p;
        tmpArray[1] = v;
        tmpArray[2] = t;
    } else if (tmp == 3) {
        tmpArray[0] = p;
        tmpArray[1] = q;
        tmpArray[2] = v;
    } else if (tmp == 4) {
        tmpArray[0] = t;
        tmpArray[1] = p;
        tmpArray[2] = v;
    } else {
        tmpArray[0] = v;
        tmpArray[1] = p;
        tmpArray[2] = q;
    }

    res[0] = (unsigned char)(tmpArray[0] * 256);
    res[1] = (unsigned char)(tmpArray[1] * 256);
    res[2] = (unsigned char)(tmpArray[2] * 256);
    res[3] = 255;

    return res;
}