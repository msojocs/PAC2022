#!/bin/bash
root_dir=$(dirname $0)
cd $root_dir
rm -rf hotspots hotspots.tar.gz
make all
vtune -collect hotspots -result-dir hotspots  ./main.exe
tar -zcf hotspots.tar.gz hotspots main.exe main.cpp main.o
rm -rf hotspots