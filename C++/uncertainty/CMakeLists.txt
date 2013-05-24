PROJECT (Riemanian_atlas_analysis)

cmake_minimum_required(VERSION 2.4)

FIND_PACKAGE(ITK)

IF(ITK_FOUND)
	INCLUDE(${USE_ITK_FILE})
ELSE(ITK_FOUND)
	MESSAGE("ITK NOT FOUND !")
ENDIF(ITK_FOUND)

FIND_PACKAGE( Boost 1.40 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )


SET(CXX_FILES
	main.cxx
	functions.cxx
)
SET(H_FILES
	functions.h
)

MESSAGE("CXX_files : ${CXX_FILES}")
MESSAGE("H_files : ${H_FILES}")

IF(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 OLD)
ENDIF(COMMAND cmake_policy)


ADD_EXECUTABLE(Riemanian_atlas_analysis ${CXX_FILES} ${H_FILES})
IF(NOT WIN32 OR CYGWIN)
  SET(BOOST_LIBRARIES boost_program_options)
ENDIF(NOT WIN32 OR CYGWIN)
TARGET_LINK_LIBRARIES(Riemanian_atlas_analysis ${ITK_LIBRARIES} ${Boost_LIBRARIES})
