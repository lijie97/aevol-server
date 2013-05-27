if WITH_IN2P3
  SUBDIRS = src
else # WITHOUT_IN2P3
  SUBDIRS = src doc
endif # IN2P3

EXTRA_DIST = src/misc

TESTS = tests/basic.in tests/space.in
TEST_EXTENSIONS = .in
IN_LOG_COMPILER = tests/simple.sh
