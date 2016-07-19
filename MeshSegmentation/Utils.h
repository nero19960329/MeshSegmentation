#pragma once

extern double distanceBetweenVectors(unsigned char a[4], unsigned char b[4]);
extern bool equals(unsigned char a[4], unsigned char b[4], double k = 1.0);
extern unsigned char* HSVtoRGB(double h, double s, double v);