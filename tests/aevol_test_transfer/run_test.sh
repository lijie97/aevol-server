#!/bin/sh

cd first_run

../../aevol_X11_debug -n 100

cd ../load_run
cp -r ../first_run/backup .
cp -r ../first_run/tree .
cp ../first_run/stat* .

../../aevol_X11_debug -f backup/gen_000050.ae -n 50

cd ..
pass=$(diff first_run/stat_fitness_glob.out load_run/stat_fitness_glob.out | wc -l)
if test "$pass" -eq "0";
  then echo "OK";
  else echo "fail";
fi
