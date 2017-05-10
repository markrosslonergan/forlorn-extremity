cmake_minimum_required (VERSION 2.6)
project (miniboone)
# The version number.
set (MINIBOONE_VERSION_MAJOR 2)
set (MINIBOONE_VERSION_MINOR 0)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
          set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
    else ()
	    set (CMAKE_CXX_STANDARD 11)
endif ()


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})



find_package(GSL REQUIRED)
if (GSL_FOUND)
	  include_directories(${GSL_INCLUDE_DIRS})
	  link_libraries (${GSL_LIBRARIES})
endif (GSL_FOUND)


# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")
add_subdirectory(src)

