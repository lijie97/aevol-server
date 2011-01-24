
../aevol_X11 -n 100 > /dev/null

../post_treatments/lineage -f backup/gen_000100.ae
../post_treatments/ancstats -f lineage-b000000-e000100-r1000.ae
../post_treatments/create_eps -f backup/gen_000100.ae

echo "OK"
  