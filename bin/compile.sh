#!/bin/bash
# Simple one command compiler
# Assumes cmake has generated an appropreate Makefile

cd ../build
make -j4 #Helpful when we have long complile chains.
make install
