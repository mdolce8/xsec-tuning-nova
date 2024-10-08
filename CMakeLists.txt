cmake_minimum_required(VERSION 3.26)
project(xsec_tuning_nova)

set(CMAKE_CXX_STANDARD 17)

include_directories(/Users/michaeldolce/Development/novasoft-code/novasoft)

# Jan. 2024: had to add this myself, bc of all the stupid Command Line Tools updates, I just set it to the most recent MacOS path.
include_directories(/Library/Developer/CommandLineTools/SDKs/MacOSX14.2.sdk/usr/include/c++/v1)

include_directories(/cvmfs/nova.opensciencegrid.org/externals/cafanacore/v01.35/include)
include_directories(/cvmfs/nova.opensciencegrid.org/externals/osclib/v00.04/include)
include_directories(/cvmfs/nova.opensciencegrid.org/externals/root/v6_16_00/Linux64bit+2.6-2.12-e17-debug/include)
#include_directories(/cvmfs/nova.opensciencegrid.org/externals/novarwgt/v00.44/inc)

include_directories(/cvmfs/nova-development.opensciencegrid.org/novasoft/releases/development/include)


FILE(GLOB_RECURSE SOURCES *.cxx *.cc)
FILE(GLOB_RECURSE HEADERS *.h)
FILE(GLOB_RECURSE MACROS *.C)


# library and target (xsec_tuning_nova) must be identical.
add_library(xsec_tuning_nova SHARED ${SOURCES} ${HEADERS} ${MACROS})
set_target_properties(xsec_tuning_nova PROPERTIES LINKER_LANGUAGE CXX)