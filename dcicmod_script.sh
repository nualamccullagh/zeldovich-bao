#!/bin/bash
gcc -I/usr/include/python2.6/ -fPIC -fopenmp -c cicdensmodule.c

gcc -shared -Xlinker -export-dynamic -I/usr/include/python2.6/ -L/usr/lib/python2.6 -lm cicdensmodule.o -o cicdens.so

