######################### SOURCE FILES ###############################
function(source_files_g2g)
    # Include source files: Common files to all kind of compilations
    file(GLOB G2G_FILES *.cpp cpu/*.cpp analytic_integral/*.cpp excited/saverho.cpp)

    # Include sorce files according to the usage or not of external libraries
    if (EXTERNAL)
        list(APPEND G2G_FILES excited/calc_FXC.cpp excited/calc_gradients.cpp
                         excited/calc_VXC.cpp excited/g2g_calcgradXC.cpp
                         excited/g2g_calculateXC.cpp excited/g2g_open_calculateXC.cpp
                         excited/recalc_densGS.cpp cioverlap/cioverlap_fake.cpp
                         libint/g2g_libint.cpp libint/libintproxy.cpp)
    else()
        list(APPEND G2G_FILES cioverlap/cioverlap_fake.cpp excited/fake_subs.cpp
                         libint/g2g_libint_fake.cpp)
    endif()

   set(SRC_FILES ${G2G_FILES} PARENT_SCOPE)
endfunction(source_files_g2g)

######################### SOURCE FILES ###############################
# openMP
function(check_openmp)
    find_package(OpenMP)
    if (NOT OPENMP_FOUND)
        message(FATAL_ERROR "OpenMP library was not found")
    else()
        target_link_libraries(g2g PUBLIC OpenMP::OpenMP_CXX)
    endif()
endfunction(check_openmp)

# External Libraries
function(link_external_libs)
    # Setting the LD_LIBRARY_PATH
    string(REPLACE ":" ";" RUNTIME_PATH "$ENV{LD_LIBRARY_PATH}")

    # LIBXC
    find_library(mylibxc NAME xc HINTS ${RUNTIME_PATH} DOC "LibXC Library")
    if (${mylibxc} STREQUAL mylibxc-NOTFOUND)
        message(FATAL_ERROR "LibXC Library was not found")
    else()
        message(STATUS "LibXC was found as ${mylibxc}")
    endif()

    # LIBINT
    find_library(mylibint NAME int2 HINTS ${RUNTIME_PATH} DOC "LibINT Library")
    if (${mylibint} STREQUAL mylibint-NOTFOUND)
        message(FATAL_ERROR "LibINT Library was not found")
    else()
        message(STATUS "LibINT was found as ${mylibint}")
    endif()

    target_compile_definitions(g2g PRIVATE "USE_LIBXC=1" "LIBXC_CPU=1" "USE_LIBINT=1")
    target_link_libraries(g2g PUBLIC ${mylibxc} ${mylibint})
    target_include_directories(g2g PRIVATE libint libxc
                                           $ENV{LIBINT_HOME}/include 
                                           $ENV{LIBINT_HOME}/include/libint2
                                           $ENV{EIGEN_HOME}
                                           $ENV{LIBXC_HOME_CPU}/include)
endfunction(link_external_libs)

function(download_external_libs)
    include(ExternalProject)
    # LIBINT: DOWNLOAD AND INSTALLATION
    set(libintHOME ${PROJECT_BINARY_DIR}/external/libintHOME)
    ExternalProject_Add(libintDownloaded
        PREFIX ${libintHOME}
        GIT_REPOSITORY https://github.com/evaleev/libint.git
        SOURCE_DIR ${libintHOME}/libint-source
        BINARY_DIR ${libintHOME}/libint-build
        CONFIGURE_COMMAND
            COMMAND cd ${libintHOME}/libint-source && ./autogen.sh
            COMMAND ${libintHOME}/libint-source/configure --prefix=${libintHOME}/libint-install --enable-eri2=1 --enable-eri=1 --enable-1body=1 --with-gnu-ld --enable-shared=yes --enable-static=no
    )
    # LIBXC: DOWNLOAD AND INSTALLATION
    set(libxcHOME ${PROJECT_BINARY_DIR}/external/libxcHOME)
    ExternalProject_Add(libxcDownloaded
        PREFIX ${libxcHOME}
        GIT_REPOSITORY https://gitlab.com/libxc/libxc.git
        GIT_TAG master
        SOURCE_DIR ${libxcHOME}/libxc-source
        BINARY_DIR ${libxcHOME}/libxc-build
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${libxcHOME}/libxc-install
                   -DBUILD_SHARED_LIBS=ON -DDISABLE_KXC=OFF -DDISABLE_LXC=OFF
    )

    # LIBINT LINKING
    set(LIBINT_INCLUDE_DIR ${libintHOME}/libint-install/include)
    set(LIBINT_LIB_FILE ${libintHOME}/libint-install/lib/${CMAKE_SHARED_LIBRARY_PREFIX}int2${CMAKE_SHARED_LIBRARY_SUFFIX})
    add_library(myLIBINT SHARED IMPORTED)
    set_target_properties(myLIBINT PROPERTIES IMPORTED_LOCATION ${LIBINT_LIB_FILE})
    include_directories(${LIBINT_INCLUDE_DIR} ${LIBINT_INCLUDE_DIR}/libint2 $ENV{EIGEN_HOME})
    target_link_libraries(g2g PUBLIC myLIBINT)
    add_dependencies(g2g libintDownloaded)

    # TODO WITH EIGEN
    #target_compile_definitions(g2g PRIVATE "USE_LIBXC=1" "LIBXC_CPU=1" "USE_LIBINT=1")
    #target_include_directories(g2g PRIVATE libint libxc
    #                                       $ENV{LIBINT_HOME}/include 
    #                                       $ENV{LIBINT_HOME}/include/libint2
    #                                       $ENV{EIGEN_HOME}
    #                                       $ENV{LIBXC_HOME_CPU}/include)

    # LIBXC LINKING
    set(LIBXC_INCLUDE_DIR ${libxcHOME}/libxc-install/include)
    set(LIBXC_LIB_FILE ${libxcHOME}/libxc-install/lib/${CMAKE_SHARED_LIBRARY_PREFIX}xc${CMAKE_SHARED_LIBRARY_SUFFIX})
    add_library(myLIBXC SHARED IMPORTED)
    set_target_properties(myLIBXC PROPERTIES IMPORTED_LOCATION ${LIBXC_LIB_FILE})
    include_directories(${LIBXC_INCLUDE_DIR})
    target_link_libraries(g2g PUBLIC myLIBXC)
    target_compile_definitions(g2g PRIVATE "USE_LIBXC=1" "LIBXC_CPU=1" "USE_LIBINT=1")
    add_dependencies(g2g libxcDownloaded)
    
endfunction(download_external_libs)
