#!/bin/sh

../aevol_X11_debug -n 100

../post_treatments/lineage_debug -f backup/gen_000100.ae
../post_treatments/ancstats_debug -f lineage-b000000-e000100.ae
../post_treatments/create_eps_debug -f backup/gen_000100.ae

echo "OK"
  