#!/bin/ksh
#
#==>  WAM post-processing pspec
#
#     Arno Behrens (May 2012)
#
#PBS -S /bin/sh
#PBS -N pspec
#PBS -j oe
#PBS -l select=1:ncpus=1
#PBS -l place=scatter
#PBS -q ksd-test
#
set -k
WAMDIR=/data/behrens/mywave
WORK=/data/behrens/mywave/workdir
cd ${WORK}/tempsg
cp ${WAMDIR}/const/Spectra_User .
cp ${WAMDIR}/abs/pspec pspec.exe
#
./pspec.exe
mv Spectra_Prot ${WAMDIR}/dayfiles/pspec_prot_c
rm Spectra_User pspec.exe
#
# ===================================================================
#  GRID FILES HAVE BEEN CREATED AND SAVED.
#  END OF JOB PGRID.
# ===================================================================
#
