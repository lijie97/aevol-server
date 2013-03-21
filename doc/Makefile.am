if HAVE_DOXYGEN

all: doxygen

doxygen: ../src/*.cpp ../src/lib-aevol/*.h ../src/lib-aevol/*.cpp ../src/post_treatments/*.cpp
	doxygen aevol.conf ; echo "document generated\n" > doxygen

latex: doxygen

html: doxygen

clean-local:
	rm -rf html latex doxygen
  
EXTRA_DIST = aevol.conf index.html html latex param.in.example

endif
