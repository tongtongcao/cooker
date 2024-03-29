# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)
# Needed for cmake 2.6
IF(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

###########################################
# POI:
# Set project name
PROJECT(cookerrawtree)

#############################################
# POI:
# Set the library version
SET(COOKERTREE_MAJOR_VERSION 1)
SET(COOKERTREE_MINOR_VERSION 1)
SET(COOKERTREE_PATCH_VERSION 0)
SET(COOKERTREE_VERSION "${COOKERTREE_MAJOR_VERSION}.${COOKERTREE_MINOR_VERSION}.${COOKERTREE_PATCH_VERSION}")
LIST(APPEND COOKERTREE_LIBRARY_PROPERTIES VERSION "${COOKERTREE_VERSION}" SOVERSION "${COOKERTREE_MAJOR_VERSION}")

# find source and header files, set paths, includes...
SET(LIBNAME "cookerrawtree")
FILE(GLOB SRCS src/*.cpp)

###################################################3
# POI:
# set header files included

FILE(GLOB HEADERS src/include/*.h
${CMAKE_SOURCE_DIR}/include/cookerrawtree.h 
${CMAKE_SOURCE_DIR}/include/mctree.h
${CMAKE_SOURCE_DIR}/include/generatortree.h
${CMAKE_SOURCE_DIR}/include/csitree.h
${CMAKE_SOURCE_DIR}/include/tof1tree.h
${CMAKE_SOURCE_DIR}/include/tof2tree.h
${CMAKE_SOURCE_DIR}/include/mwpctree.h
${CMAKE_SOURCE_DIR}/include/actree.h
${CMAKE_SOURCE_DIR}/include/anaE36.h
${CMAKE_SOURCE_DIR}/include/cookedTof.h
${CMAKE_SOURCE_DIR}/include/trackingE36.h
${CMAKE_SOURCE_DIR}/include/pgctree.h
${CMAKE_SOURCE_DIR}/include/mergedCaliTargetTofMwpc.h
${CMAKE_SOURCE_DIR}/include/track.h
${CMAKE_SOURCE_DIR}/include/trackArray.h
)

SET(INCLUDES ${ROOT_INCLUDE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/src/include" "${CMAKE_CURRENT_SOURCE_DIR}/../../include")
INCLUDE_DIRECTORIES(${INCLUDES})

# create root dictionary
SET(LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/src/include/LinkDef.hh)
SET(DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/${LIBNAME}Dict.cpp)
ROOT_GENERATE_DICTIONARY("${HEADERS}" "${LINKDEF}" "${DICTIONARY}" "${INCLUDES}")
LIST(APPEND SRCS ${DICTIONARY})

# build the library
ADD_LIBRARY(${LIBNAME} SHARED ${SRCS})
SET(ROOTSTUFF Gui Gpad Hist Graf Graf3d Tree Rint Postscript Matrix Physics MathCore RIO Net Thread Core Cint Eve Ged RGL TreePlayer EG Geom)
TARGET_LINK_LIBRARIES(${LIBNAME}  "-L${ROOT_LIBRARY_DIR}" ${ROOTSTUFF})

SET_TARGET_PROPERTIES(${LIBNAME} PROPERTIES ${COOKERTREE_LIBRARY_PROPERTIES})


# install the library
INSTALL(TARGETS ${LIBNAME} RUNTIME DESTINATION ${BINDIR} LIBRARY DESTINATION ${LIBDIR} ARCHIVE DESTINATION ${LIBDIR})

IF(INSTALL_HEADER)
	INSTALL(FILES ${HEADERS} DESTINATION ${INCDIR})
ENDIF(INSTALL_HEADER)

