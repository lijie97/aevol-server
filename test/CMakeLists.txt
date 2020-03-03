# Google Test headers
add_subdirectory(googletest)
include_directories("${gtest_SOURCE_DIR}/include")

# We need zlib
find_package(ZLIB REQUIRED)

# Libraries to link unit tests against
set(test_libs
  aevol
  gtest_main
)

# List of all the tests
set(TESTS
  SampleTest
  IndividualTest
  JumpingMTTest
  TranscriptionTranslationTest
)

# Create a runner for each unit test. Generates a list named TEST_RUNNERS
foreach (TEST IN LISTS TESTS)

  set(TEST_RUNNER run_${TEST})
  set(TEST_RUNNERS ${TEST_RUNNERS} ${TEST_RUNNER})

  add_executable(${TEST} ${TEST}.cpp)

  set_target_properties(${TEST}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
  )

  target_link_libraries(${TEST} ${test_libs})

  add_custom_target(${TEST_RUNNER}
    COMMENT "Running test ${TEST}..."
    COMMAND ${CMAKE_BINARY_DIR}/test/${TEST}
    DEPENDS ${TEST}
  )

  set_target_properties(${TEST_RUNNER}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
  )
endforeach(TEST)

# Copy data files for TranscriptionTranslationTest
add_custom_command(TARGET TranscriptionTranslationTest
  PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/test/TranscriptionTranslationTest_files
    ${CMAKE_BINARY_DIR}/test/TranscriptionTranslationTest_files
)


# Add target check (main target for tests)
add_custom_target(check DEPENDS ${TEST_RUNNERS})
