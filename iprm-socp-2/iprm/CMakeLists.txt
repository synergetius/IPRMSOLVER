cmake_minimum_required(VERSION 3.18)

#
#set(IPRM_BUILD_TYPE Debug)

set(IPRM_VERSION_MAJOR "0")
set(IPRM_VERSION_MINOR "1")
set(IPRM_VERSION_PATCH "4")
set(IPRM_VERSION "${IPRM_VERSION_MAJOR}.${IPRM_VERSION_MINOR}.${IPRM_VERSION_PATCH}")

# Project name
project(iprm VERSION ${IPRM_VERSION})
message(STATUS "Building IPRM v${IPRM_VERSION}")

# If build type is not specified set to Release.
if(IPRM_BUILD_TYPE STREQUAL Debug)
    set(IPRM_BUILD_TYPE ${IPRM_BUILD_TYPE})
    add_compile_definitions(IPRM_DEBUG)
    set(CMAKE_C_FLAGS "-g -Wall")
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux" OR ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address,undefined")
    endif()
else()
    set(IPRM_BUILD_TYPE Release)
    set(CMAKE_C_FLAGS "-O3 -Wall")
endif()

# Detect operating system.
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    add_compile_definitions(IS_LINUX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra")
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    add_compile_definitions(IS_MACOS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra")
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    add_compile_definitions(IS_WINDOWS)
    if(IPRM_BUILD_TYPE STREQUAL Release)
        #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Ox")
        #此处替换一下
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3") 
    endif()
endif()

message(STATUS "Build Type: " ${IPRM_BUILD_TYPE})
message(STATUS "Build Flags: " ${CMAKE_C_FLAGS})

# Set integer size to 32 bits.
message(STATUS "Using 32 byte integers")
set(QDLDL_LONG OFF)

# Set floating points to double precision.
message(STATUS "Using double precision floating point")


# Add -fPIC option explicity.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add qdldl.
add_subdirectory(lib/qdldl)

# Add amd
add_subdirectory(lib/amd)

# Are we building for matlab. If so we will define printf as mexPrintf to print into matlab terminal.
if(${MATLAB})
    add_compile_definitions(MATLAB)
endif()

# Build iprm shared library.
add_library(iprm SHARED)
target_link_libraries(iprm qdldl amd)
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux" OR ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    target_link_libraries(iprm m)
endif()
target_include_directories(iprm PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(iprm PUBLIC lib/amd)
target_sources(iprm PRIVATE src/iprm_api.c src/input_validation.c src/linalg.c src/kkt.c src/utils.c src/cone.c src/iprm_error.c )
target_compile_definitions(iprm PUBLIC "IPRM_VERSION_MAJOR=${IPRM_VERSION_MAJOR}")
target_compile_definitions(iprm PUBLIC "IPRM_VERSION_MINOR=${IPRM_VERSION_MINOR}")
target_compile_definitions(iprm PUBLIC "IPRM_VERSION_PATCH=${IPRM_VERSION_PATCH}")

# Build iprm static library.
add_library(iprmstatic STATIC)
target_link_libraries(iprmstatic qdldlstatic amd)
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux" OR ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    target_link_libraries(iprmstatic m)
endif()
target_include_directories(iprmstatic PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(iprmstatic PUBLIC lib/amd)
target_sources(iprmstatic PRIVATE src/iprm_api.c src/input_validation.c src/linalg.c src/kkt.c src/utils.c src/cone.c src/iprm_error.c )
target_compile_definitions(iprmstatic PUBLIC "IPRM_VERSION_MAJOR=${IPRM_VERSION_MAJOR}")
target_compile_definitions(iprmstatic PUBLIC "IPRM_VERSION_MINOR=${IPRM_VERSION_MINOR}")
target_compile_definitions(iprmstatic PUBLIC "IPRM_VERSION_PATCH=${IPRM_VERSION_PATCH}")

# Include timer file depending on OS.
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    target_sources(iprm PRIVATE src/timer_linux.c)
    target_sources(iprmstatic PRIVATE src/timer_linux.c)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    target_sources(iprm PRIVATE src/timer_macos.c)
    target_sources(iprmstatic PRIVATE src/timer_macos.c)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    target_sources(iprm PRIVATE src/timer_windows.c)
    target_sources(iprmstatic PRIVATE src/timer_windows.c)
endif()



# Build iprm demo
#add_executable(iprm_demo examples/iprm_demo.c)
#target_include_directories(iprm_demo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_link_libraries(iprm_demo iprmstatic)

add_executable(demo_socp examples/demo_socp.c)
target_include_directories(demo_socp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(demo_socp iprmstatic)

# Test definitions.
if(ENABLE_TESTING)
    configure_file(CMakeLists.txt.in
            googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
    add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
            ${CMAKE_BINARY_DIR}/googletest-build)
    enable_testing()
    add_subdirectory(tests)
endif()

install(
    TARGETS iprmstatic
    EXPORT iprmstatic
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
install(
    TARGETS iprm
    EXPORT iprm
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
install(FILES ${iprm_headers}
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/iprm"
)
