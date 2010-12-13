#!/bin/sh

cd first_run

if test "$1" = "x"
then
  ../../aevol_X11 -n 100 > /dev/null
else
  if test "$1" = "nox"
  then
    ../../aevol -n 100 > /dev/null
  else
    echo "  usage : ./run_test x|nox"
    exit 1
  fi
fi

cd ../load_run
cp -r ../first_run/backup .
cp -r ../first_run/tree .
cp ../first_run/stat* .

if test "$1" = "x"
then
  ../../aevol_X11 -f backup/gen_000050.ae -n 50 > /dev/null
else
  ../../aevol -f backup/gen_000050.ae -n 50 > /dev/null
fi

cd ..
pass=$(diff first_run/stat_fitness_glob.out load_run/stat_fitness_glob.out | wc -l)
if test "$pass" -eq "0";
  then echo "OK";
  else echo "fail";
fi
