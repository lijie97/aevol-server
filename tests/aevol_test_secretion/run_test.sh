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

echo "OK"
