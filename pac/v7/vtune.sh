#!/bin/bash
root_dir=$(dirname $0)
version=$(basename $root_dir)

type="hotspots"
# type="uarch-exploration"
# type="memory-access"

cd $root_dir
rm -rf "${type}_${version}" "${type}_${version}.tar.gz"
make all
vtune -collect $type -result-dir "${type}_${version}"  ./main.exe
tar -zcf "${type}_${version}.tar.gz" "${type}_${version}" main.exe main.cpp main.o
rm -rf "${type}_${version}"