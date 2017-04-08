#!/bin/sh
#g++ -O4 -g svdDynamic.c RayTracer.c utils.c -lm -o RayTracer
g++ -g -O4 -fopenmp -std=gnu++11 svdDynamic.c RayTracer.c utils.c -lm -o RayTracer

