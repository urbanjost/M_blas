#!/bin/bash
export FPM_COMPILER
for PROFILE in debug release
do
# take out or add the compilers you have
for FPM_COMPILER in gfortran ifort nvfortran 
do
   (
   exec 2>&1
   time (
      exec 2>&1
      fpm  test  cblat1    -profile  $PROFILE
      fpm  test  cblat2    -profile  $PROFILE  <  test/cblat2.in
      fpm  test  cblat2_a  -profile  $PROFILE  <  test/cblat2.in
      fpm  test  cblat3    -profile  $PROFILE  <  test/cblat3.in
      fpm  test  dblat1    -profile  $PROFILE
      fpm  test  dblat1_a  -profile  $PROFILE
      fpm  test  dblat2    -profile  $PROFILE  <  test/dblat2.in
      fpm  test  dblat3    -profile  $PROFILE  <  test/dblat3.in
      fpm  test  dblat3_a  -profile  $PROFILE  <  test/dblat3.in
      fpm  test  sblat1    -profile  $PROFILE
      fpm  test  sblat2    -profile  $PROFILE  <  test/sblat2.in
      fpm  test  sblat3    -profile  $PROFILE  <  test/sblat3.in
      fpm  test  zblat1    -profile  $PROFILE
      fpm  test  zblat2    -profile  $PROFILE  <  test/zblat2.in
      fpm  test  zblat3    -profile  $PROFILE  <  test/zblat3.in
      
      cat cblat2.out
      cat cblat3.out
      cat dblat2.out
      cat dblat3.out
      cat sblat2.out
      cat sblat3.out
      cat zblat2.out
      cat zblat3.out
      
      rm -f cblat2.out cblat3.out dblat2.out dblat3.out sblat2.out sblat3.out zblat2.out zblat3.out
   
   )
   )
   done|cat -n|tee run.log.$FPM_COMPILER
   
   grep -i fail run.log.$FPM_COMPILER
done
exit
