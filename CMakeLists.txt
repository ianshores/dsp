cmake_minimum_required(VERSION 3.16)
project(ShoresDSP VERSION 0.2.0 LANGUAGES CXX)

add_library(shores_dsp INTERFACE)
add_library(shores_dsp::shores_dsp ALIAS shores_dsp)
target_include_directories(shores_dsp INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(shores_dsp INTERFACE cxx_std_17)

include(CTest)
if(BUILD_TESTING)
  add_executable(dsp_biquad_test tests/biquad.cpp)
  target_link_libraries(dsp_biquad_test PRIVATE shores_dsp)
  add_test(NAME biquad COMMAND dsp_biquad_test)
  foreach(mode native scalar sse half packed)
    add_executable(dsp_fft_${mode} tests/fft.cpp)
    target_link_libraries(dsp_fft_${mode} PRIVATE shores_dsp)
    if(mode STREQUAL "scalar")
      target_compile_definitions(dsp_fft_${mode} PRIVATE DSP_FFT_SCALAR)
    elseif(mode STREQUAL "sse")
      target_compile_definitions(dsp_fft_${mode} PRIVATE DSP_DISABLE_AVX)
    elseif(mode STREQUAL "half")
      target_compile_definitions(dsp_fft_${mode} PRIVATE DSP_TEST_HALF_FORWARD)
    elseif(mode STREQUAL "packed")
      target_compile_definitions(dsp_fft_${mode} PRIVATE DSP_TEST_PACKED DSP_TEST_PACKED8)
    endif()
    add_test(NAME fft_${mode} COMMAND dsp_fft_${mode})
  endforeach()
endif()
option(DSP_BUILD_BENCHMARKS "Build standalone DSP benchmarks" OFF)
if(DSP_BUILD_BENCHMARKS)
  add_executable(dsp_biquad_benchmark tests/biquad_benchmark.cpp)
  target_link_libraries(dsp_biquad_benchmark PRIVATE shores_dsp)
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
install(TARGETS shores_dsp EXPORT ShoresDSPTargets)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT ShoresDSPTargets NAMESPACE shores_dsp::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ShoresDSP)
configure_package_config_file(cmake/ShoresDSPConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/ShoresDSPConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ShoresDSP)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/ShoresDSPConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ShoresDSPConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/ShoresDSPConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ShoresDSP)
