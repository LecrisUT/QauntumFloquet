#!/bin/bash
echo $DYLD_LIBRARY_PATH
echo $MKLROOT
source /opt/intel/oneapi/setvars.sh --force
echo $DYLD_LIBRARY_PATH
echo $MKLROOT
ps -p $$
echo $SHELL