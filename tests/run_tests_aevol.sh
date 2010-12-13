#!/bin/sh

cd ..
./configure
make
cd tests

echo -n "running aevol_test_align... "
cd ./aevol_test_align
./run_test.sh x
cd ..

echo -n "running aevol_test_backup... "
cd ./aevol_test_backup
./run_test.sh x
cd ..

echo -n "running aevol_test_env-var... "
cd ./aevol_test_env-var
./run_test.sh x
cd ..

echo -n "running aevol_test_neutral-zones... "
cd ./aevol_test_neutral-zones
./run_test.sh x
cd ..

echo -n "running aevol_test_plasmids... "
cd ./aevol_test_plasmids
./run_test.sh x
cd ..

echo -n "running aevol_test_secretion... "
cd ./aevol_test_secretion
./run_test.sh x
cd ..