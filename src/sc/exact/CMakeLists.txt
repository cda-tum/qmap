#
# This file is part of the MQT QMAP library released under the MIT license. See README.md or go to
# https://github.com/cda-tum/mqt-qmap for more information.
#

set(MQT_QMAP_SC_EXACT_TARGET_NAME ${MQT_QMAP_TARGET_NAME}-sc-exact)

if(NOT TARGET ${MQT_QMAP_SC_EXACT_TARGET_NAME})
  # collect headers and source files
  file(GLOB SC_EXACT_HEADERS ${MQT_QMAP_INCLUDE_BUILD_DIR}/sc/exact/*.hpp)
  file(GLOB SC_EXACT_SOURCES **.cpp)

  # add SC Exact Package library
  add_library(${MQT_QMAP_SC_EXACT_TARGET_NAME} ${SC_EXACT_HEADERS} ${SC_EXACT_SOURCES})

  # set include directories
  target_include_directories(${MQT_QMAP_SC_EXACT_TARGET_NAME}
                             PUBLIC $<BUILD_INTERFACE:${MQT_QMAP_INCLUDE_BUILD_DIR}>)

  # link to the MQT::Core libraries
  target_link_libraries(
    ${MQT_QMAP_SC_EXACT_TARGET_NAME}
    PUBLIC MQT::QMapSC
    PRIVATE MQT::LogicBlocks MQT::ProjectWarnings MQT::ProjectOptions)

  # add MQT alias
  add_library(MQT::QMapSCExact ALIAS ${MQT_QMAP_SC_EXACT_TARGET_NAME})
endif()
