#!/bin/bash

gcc -dynamiclib -I/anaconda/include/python2.7/ -L/anaconda/lib/python2.7/ -lm -lpython2.7 -o cicdens.dylib cicdensmodule.c
mv cicdens.dylib cicdens.so