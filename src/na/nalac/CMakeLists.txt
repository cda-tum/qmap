set(MQT_QMAP_NALAC_TARGET_NAME ${MQT_QMAP_TARGET_NAME}-nalac)

if(NOT TARGET ${MQT_QMAP_NALAC_TARGET_NAME})
  file(GLOB NALAC_HEADERS ${MQT_QMAP_INCLUDE_BUILD_DIR}/na/nalac/*.hpp)
  file(GLOB NALAC_SOURCES *.cpp)

  add_library(${MQT_QMAP_NALAC_TARGET_NAME} ${NALAC_HEADERS} ${NALAC_SOURCES})

  target_include_directories(${MQT_QMAP_NALAC_TARGET_NAME}
                             PUBLIC $<BUILD_INTERFACE:${MQT_QMAP_INCLUDE_BUILD_DIR}>)

  target_link_libraries(
    ${MQT_QMAP_NALAC_TARGET_NAME}
    PUBLIC MQT::QMapNA
    PRIVATE MQT::ProjectOptions MQT::ProjectWarnings)

  add_library(MQT::NALAC ALIAS ${MQT_QMAP_NALAC_TARGET_NAME})
endif()
