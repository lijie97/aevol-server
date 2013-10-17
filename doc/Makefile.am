
all:

if HAVE_DOXYGEN

doxygen: ../src/*.cpp ../src/libaevol/*.h ../src/libaevol/*.cpp ../src/post_treatments/*.cpp
	doxygen aevol.conf ; echo "document generated\n" > doxygen

latex: doxygen

html: doxygen

doc: doxygen

clean: clean-local

clean-local:
	rm -rf html doxygen

endif

dist_man1_MANS = aevol_create.1 aevol_run.1 aevol_run_X11.1 aevol_modify.1 aevol_propagate.1
EXTRA_DIST = README.txt aevol.conf
