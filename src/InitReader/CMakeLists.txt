# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)
# Needed for cmake 2.6
IF(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

# Set project name
PROJECT(InitReader)

SET(INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/include" ${Boost_INCLUDE_DIRS} "${CMAKE_SOURCE_DIR}/include" "${XERCESC_INCLUDE_DIR}")

INCLUDE_DIRECTORIES(${INCLUDES})

add_library(InitReader src/InitReader.cpp )
TARGET_LINK_LIBRARIES(InitReader xqilla xerces-c xerces-c-3.1)
