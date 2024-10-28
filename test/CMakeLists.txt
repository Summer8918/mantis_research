include_directories("include;/home/jie/library/sdsl-lite/sdsl_install_dir/include;/home/jie/library/sdsl-lite/sdsl_install_dir/lib")
link_directories("/home/jie/library/sdsl-lite/sdsl_install_dir/include;/home/jie/library/sdsl-lite/sdsl_install_dir/lib")
include_directories("../include")
include_directories("../src")

add_executable(test_hash test_hash.cc ../src/gqf/hashutil.c)
add_executable(sdsl_test sdsl_test.cc)
add_executable(validateMantisUtilsTest validateMantisUtilsTest.cc ../src/validateMantisUtils.cc)

target_link_libraries(test_hash sdsl Threads::Threads)
target_link_libraries(sdsl_test sdsl)

if (SDSL_INSTALL_PATH)
   set_property(TARGET test_hash APPEND_STRING PROPERTY LINK_FLAGS "-L${SDSL_INSTALL_PATH}/lib")
endif()


