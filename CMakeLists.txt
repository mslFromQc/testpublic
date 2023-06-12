cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS /W3 /GR /EHsc /DUNICODE /Zc:__cplusplus")
set(CMAKE_CXX_FLAGS_DEBUG "/MDd /Zi /Ob0 /Od /RTC1")
set(CMAKE_CXX_FLAGS_RELEASE "/MD /O2 /Ob2 /DNDEBUG")

add_compile_definitions(RAPIDJSON_HAS_STDSTRING=1)

include_directories(.)
include_directories(components/hdf5/include)
include_directories(components/rapidjson/include)
 
link_directories(components/hdf5/lib/${CMAKE_BUILD_TYPE})

file(GLOB
  _FileFormatChanges
  ${CMAKE_CURRENT_SOURCE_DIR}/compatibilities/*.cpp
)

add_executable("create-file" examples/create-file/main.cpp)
target_link_libraries("create-file" optimized libhdf5 optimized zlibstatic debug libhdf5_D debug zlibstaticd)

add_executable("read-setup" examples/read-setup/main.cpp)
target_link_libraries("read-setup" optimized libhdf5 optimized zlibstatic debug libhdf5_D debug zlibstaticd)

add_executable("write-ascan" examples/write-ascan/main.cpp)
target_link_libraries("write-ascan" optimized libhdf5 optimized zlibstatic debug libhdf5_D debug zlibstaticd)

add_executable("write-cscan" examples/write-cscan/main.cpp)
target_link_libraries("write-cscan" optimized libhdf5 optimized zlibstatic debug libhdf5_D debug zlibstaticd)

add_executable("read-ascan" examples/read-ascan/main.cpp ${_FileFormatChanges})
target_link_libraries("read-ascan" optimized libhdf5 optimized zlibstatic debug libhdf5_D debug zlibstaticd)

add_executable("version-compatibilities" examples/version-compatibilities/main.cpp ${_FileFormatChanges})
target_link_libraries("version-compatibilities" optimized libhdf5 optimized zlibstatic debug libhdf5_D debug zlibstaticd)

