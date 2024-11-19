#!/bin/sh 

echo conda info -s

conda list | grep json-fortran
which gfortran