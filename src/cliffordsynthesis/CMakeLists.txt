#
# This file is part of the MQT QMAP library released under the MIT license. See README.md or go to
# https://github.com/cda-tum/mqt-qmap for more information.
#

set(MQT_QMAP_CLIFFORD_SYNTHESIS_TARGET_NAME ${MQT_QMAP_TARGET_NAME}-cliffordsynthesis)

if(NOT TARGET ${MQT_QMAP_CLIFFORD_SYNTHESIS_TARGET_NAME})
  # collect headers and source files
  file(GLOB_RECURSE CLIFFORD_SYNTHESIS_HEADERS
       ${MQT_QMAP_INCLUDE_BUILD_DIR}/cliffordsynthesis/*.hpp)
  file(GLOB_RECURSE CLIFFORD_SYNTHESIS_SOURCES **.cpp)
  list(APPEND CLIFFORD_SYNTHESIS_SOURCES "${CMAKE_SOURCE_DIR}/src/sc/utils.cpp")

  # add CliffordSynthesis Package library
  add_library(${MQT_QMAP_CLIFFORD_SYNTHESIS_TARGET_NAME} ${CLIFFORD_SYNTHESIS_HEADERS}
                                                         ${CLIFFORD_SYNTHESIS_SOURCES})
  message(${CLIFFORD_SYNTHESIS_SOURCES})
  # set include directories
  target_include_directories(${MQT_QMAP_CLIFFORD_SYNTHESIS_TARGET_NAME}
                             PUBLIC $<BUILD_INTERFACE:${MQT_QMAP_INCLUDE_BUILD_DIR}>)

  # link to the MQT::Core libraries
  target_link_libraries(
    ${MQT_QMAP_CLIFFORD_SYNTHESIS_TARGET_NAME}
    PUBLIC MQT::CoreIR plog::plog nlohmann_json::nlohmann_json
    PRIVATE MQT::LogicBlocks MQT::ProjectWarnings MQT::ProjectOptions MQT::CoreCircuitOptimizer)

  # add MQT alias
  add_library(MQT::QMapCliffordSynthesis ALIAS ${MQT_QMAP_CLIFFORD_SYNTHESIS_TARGET_NAME})
endif()
