#[[---------------------------------------------------------------------------]
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
[-----------------------------------------------------------------------------]
License
    This file is part of a specialisation of the Approximate Nearest Neighbor
    Library (ANN) library for HELYXcore.  This software is provided under the
    provisions of the Lesser GNU Public License (LGPL), in accordance with the
    requirements of the ANN library.  See the file ../ReadMe.txt for further
    information.

Copyright
    (c) 2020 Engys Ltd.


[----------------------------------------------------------------------------]]

set(TARGET_NAME "knn")

add_helyx_library(${TARGET_NAME} SHARED
./src/ANN.C
./src/bd_fix_rad_search.C
./src/bd_pr_search.C
./src/bd_search.C
./src/bd_tree.C
./src/brute.C
./src/kd_dump.C
./src/kd_fix_rad_search.C
./src/kd_pr_search.C
./src/kd_search.C
./src/kd_split.C
./src/kd_tree.C
./src/kd_util.C
./src/perf.C
./src/knn.C
)

if("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
    target_compile_definitions(knn
        #PRIVATE WIN32
        PRIVATE DLL_EXPORTS
    )
endif()
