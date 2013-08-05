if HAVE_DOXYGEN

all: doxygen

doxygen: ../src/*.cpp ../src/libaevol/*.h ../src/libaevol/*.cpp ../src/post_treatments/*.cpp
	doxygen aevol.conf ; echo "document generated\n" > doxygen

latex: doxygen

html: doxygen

clean: clean-local

clean-local:
	rm -rf html doxygen
  
EXTRA_DIST = param_example.in README.txt aevol.conf

endif
