#!/bin/bash

DEST=tests/`basename ${1} .in`
PARAM=${1}
mkdir -p ${DEST}
cp ${PARAM} ${DEST}/param.in
cd ${DEST}

# Chech that aevol has been compiled correctly
if [ -x "../../src/aevol_run" ]; then
  RUN="../../src/aevol_run"
elif [ -x "../../src/aevol_run_debug" ]; then
  RUN="../../src/aevol_run_debug"
else
  echo "aevol_run not found"
  exit -1
fi

if [ ! -x "../../src/aevol_create" ]; then
  echo "aevol_create not found"
  exit -1
fi

# Create a simple experiment and run it for 10 generations
../../src/aevol_create -f param.in
${RUN} -r 0 -n 10

# Check that last generation has been saved
if [ ! -f "populations/pop_000010.ae" ]; then
  echo "Backup not found for generation 10"
  exit -1
fi

# Check the content of a stat file
if [ `grep '^[1-9]' stats/stat_fitness_glob.out |wc -l` != 10 ]; then
  echo "stat_fitness_glob.out does not contain the expected number of generation"
  exit -1
fi

exit 0
