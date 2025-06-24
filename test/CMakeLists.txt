# test/CMakeLists.txt

#add_subdirectory(base_lio)

if(USE_LIBINT)
    add_subdirectory(external_functionals)
    if(USE_LIBXC_CPU OR USE_LIBXC_GPU)
        add_subdirectory(linear_response)
    endif()
endif()