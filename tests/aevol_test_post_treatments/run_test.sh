if test "$1" = "x"
then
  ../aevol_X11 -n 100 > /dev/null
else
  if test "$1" = "nox"
  then
    ../aevol -n 100 > /dev/null
  else
    echo "  usage : ./run_test x|nox"
  fi
fi

../post_treatments/lineage -f backup/gen_000100.ae
../post_treatments/ancstats -f lineage-b000000-e000100-r1000.ae
../post_treatments/create_eps -f backup/gen_000100.ae

echo "OK"
  