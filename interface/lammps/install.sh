#!/bin/sh
LAMMPS_PATH=$1
curdir=$(pwd)
cp ../../src/nep.*  USER-NEP/
cp -r USER-NEP/ $LAMMPS_PATH/src
cd $LAMMPS_PATH/src
make clean-serial
make no-user-nep
make yes-user-nep
make serial
mv lmp_serial $curdir
