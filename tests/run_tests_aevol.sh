#!/bin/sh

cd ..
./configure --with-debug
make clean
make
cd tests

echo -n "running aevol_test_align... "
cd ./aevol_test_align
./run_test.sh
cd ..

echo -n "running aevol_test_transfer... "
cd ./aevol_test_transfer
./run_test.sh
cd ..

echo -n "running aevol_test_post_treatments... "
cd ./aevol_test_post_treatments
./run_test.sh
cd ..

echo -n "running aevol_test_env-var... "
cd ./aevol_test_env-var
./run_test.sh
cd ..

echo -n "running aevol_test_plasmids... "
cd ./aevol_test_plasmids
./run_test.sh
cd ..

echo -n "running aevol_test_secretion... "
cd ./aevol_test_secretion
./run_test.sh
cd ..
