if HAVE_DOXYGEN

all: doxygen

doxygen: ../src/*.cpp ../src/libaevol/*.h ../src/libaevol/*.cpp ../src/post_treatments/*.cpp
	doxygen aevol.conf ; echo "document generated\n" > doxygen

latex: doxygen

html: doxygen

clean-local:
	rm -rf html doxygen
  
EXTRA_DIST = param.in.example README.txt

endif