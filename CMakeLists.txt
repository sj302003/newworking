#==============================================================================
# NLOPT CMake file
#
# NLopt is a free/open-source library for nonlinear optimization, providing
# a common interface for a number of different free optimization routines
# available online as well as original implementations of various other
# algorithms
# WEBSITE: http://ab-initio.mit.edu/wiki/index.php/NLopt
# AUTHOR: Steven G. Johnson
#
# This CMakeLists.txt file was created to compile NLOPT with the CMAKE utility.
# Benoit Scherrer, 2010 CRL, Harvard Medical School
# Copyright (c) 2008-2009 Children's Hospital Boston
#==============================================================================
cmake_minimum_required (VERSION 3.2)

if (NOT DEFINED CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif ()

project (nlopt)

#==============================================================================
# version
set (NLOPT_MAJOR_VERSION "2")
set (NLOPT_MINOR_VERSION "7")
set (NLOPT_BUGFIX_VERSION "1")
set (NLOPT_VERSION_STRING ${NLOPT_MAJOR_VERSION}.${NLOPT_MINOR_VERSION}.${NLOPT_BUGFIX_VERSION})
message (STATUS "NLopt version ${NLOPT_VERSION_STRING}")

# This is the ABI version number, which differes from the API version above
# (it indicates ABI compatibility), but they are typically incremented together.
set(SO_MAJOR 0)
set(SO_MINOR 11)
set(SO_PATCH 1)
#==============================================================================

list (APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

option (NLOPT_CXX "enable cxx routines" ON)
option (NLOPT_FORTRAN "enable fortran tests" OFF)
option (BUILD_SHARED_LIBS "Build NLopt as a shared library" ON)
option (NLOPT_PYTHON "build python bindings" ON)
option (NLOPT_OCTAVE "build octave bindings" ON)
option (NLOPT_MATLAB "build matlab bindings" ON)
option (NLOPT_GUILE "build guile bindings" ON)
option (NLOPT_SWIG "use SWIG to build bindings" ON)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  option (NLOPT_TESTS "build unit tests" ON)
else ()
  option (NLOPT_TESTS "build unit tests" OFF)
endif ()

if (NLOPT_FORTRAN)
  enable_language (Fortran)
endif ()

include (GNUInstallDirs)

# Offer the user the choice of overriding the installation directories
set (INSTALL_LIB_DIR     ${CMAKE_INSTALL_LIBDIR}        CACHE PATH "Installation directory for libraries")
set (INSTALL_BIN_DIR     ${CMAKE_INSTALL_BINDIR}        CACHE PATH "Installation directory for executables")
set (INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}    CACHE PATH "Installation directory for header files")
set (INSTALL_DATA_DIR    ${CMAKE_INSTALL_DATADIR}/nlopt CACHE PATH "Installation directory for data files")
set (INSTALL_MAN_DIR     ${CMAKE_INSTALL_MANDIR}        CACHE PATH "Installation directory for man documentation")
set (INSTALL_CMAKE_DIR   ${INSTALL_LIB_DIR}/cmake/nlopt CACHE PATH "Installation directory for cmake config files")

# Make relative paths absolute (needed later on)
foreach (p LIB BIN INCLUDE DATA CMAKE)
  set (var INSTALL_${p}_DIR)
  set (RELATIVE_INSTALL_${p}_DIR ${INSTALL_${p}_DIR})
  if (NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif ()
endforeach ()


set (CMAKE_INSTALL_RPATH ${INSTALL_LIB_DIR})
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(POLICY CMP0042)
  # Set MACOSX_RPATH to ON
  cmake_policy(SET CMP0042 NEW)
endif()
if (POLICY CMP0086)
  # UseSWIG honors SWIG_MODULE_NAME via -module flag
  cmake_policy(SET CMP0086 NEW)
endif ()

include (CheckIncludeFiles)
include (CheckFunctionExists)
include (CheckTypeSize)
include (CheckCCompilerFlag)
include (CheckCXXSymbolExists)
include (CheckCXXCompilerFlag)
include (CheckLibraryExists)

#==============================================================================
# COMPILATION CHECKINGS and CONFIGURATION GENERATION
#==============================================================================
check_include_file (getopt.h HAVE_GETOPT_H)
check_include_file (unistd.h HAVE_UNISTD_H)
check_include_file (stdint.h HAVE_STDINT_H)
check_include_file (time.h HAVE_TIME_H)
check_include_file (sys/time.h HAVE_SYS_TIME_H)
if (HAVE_TIME_H AND HAVE_SYS_TIME_H)
  set (TIME_WITH_SYS_TIME TRUE)
endif ()
check_function_exists (getpid HAVE_GETPID)
check_function_exists (syscall HAVE_GETTID_SYSCALL)
check_function_exists (isinf HAVE_ISINF)
check_function_exists (isnan HAVE_ISNAN)
check_function_exists (gettimeofday HAVE_GETTIMEOFDAY)
check_function_exists (qsort_r HAVE_QSORT_R)
check_function_exists (time HAVE_TIME)
check_function_exists (copysign HAVE_COPYSIGN)
check_function_exists (getopt HAVE_GETOPT)
check_type_size ("uint32_t" SIZEOF_UINT32_T)
set (HAVE_UINT32_T ${SIZEOF_UINT32_T})
check_type_size ("unsigned int" SIZEOF_UNSIGNED_INT)
check_type_size ("unsigned long" SIZEOF_UNSIGNED_LONG)

check_library_exists ("m" sqrt "" HAVE_LIBM)
if (HAVE_LIBM)
  set (M_LIBRARY m)
endif()

if (NOT DEFINED HAVE_FPCLASSIFY)
  message(STATUS "Looking for fpclassify")
  file (WRITE ${PROJECT_BINARY_DIR}/fpclassify.c "#include <math.h>\n")
  file (APPEND ${PROJECT_BINARY_DIR}/fpclassify.c "int main(void) {\n")
  file (APPEND ${PROJECT_BINARY_DIR}/fpclassify.c "if (!fpclassify(3.14159)) fpclassify(2.7183);\n")
  file (APPEND ${PROJECT_BINARY_DIR}/fpclassify.c "  return 0; }\n")
  try_compile (HAVE_FPCLASSIFY
  ${PROJECT_BINARY_DIR}/build_fpclassify
  ${PROJECT_BINARY_DIR}/fpclassify.c
  CMAKE_FLAGS -DLINK_LIBRARIES=m)
  message(STATUS "Looking for fpclassify - ${HAVE_FPCLASSIFY}")
endif ()

option (WITH_THREADLOCAL "check thread local keyword" ON)
if (WITH_THREADLOCAL AND NOT DEFINED THREADLOCAL)
  foreach (_THREADLOCAL_KEY "__thread" "__declspec(thread)")
    unset (HAVE_THREAD_LOCAL_STORAGE CACHE)
    check_c_source_compiles("
    ${_THREADLOCAL_KEY} int tls;

    int main(void) {
        return 0;
    }" HAVE_THREAD_LOCAL_STORAGE)
    if (HAVE_THREAD_LOCAL_STORAGE)
      set (THREADLOCAL ${_THREADLOCAL_KEY} CACHE STRING "Thread local keyword")
    endif ()
  endforeach()
endif ()


if (NLOPT_CXX OR NLOPT_PYTHON OR NLOPT_GUILE OR NLOPT_OCTAVE)
  check_cxx_symbol_exists (__cplusplus ciso646 SYSTEM_HAS_CXX)
  if (SYSTEM_HAS_CXX)
    check_cxx_compiler_flag ("-std=c++17" SUPPORTS_STDCXX17)
    if (SUPPORTS_STDCXX17)
      set (CMAKE_CXX_FLAGS "-std=c++17 -O3 -march=native -ffast-math -fpermissive -fopenmp ${CMAKE_CXX_FLAGS}")
     # set (CMAKE_CXX_FLAGS "-std=c++17 -march=native -ffast-math -fpermissive -fopenmp ${CMAKE_CXX_FLAGS}")
      if (NLOPT_CXX)
        set (NLOPT_CXX11 ON)
      endif ()
    endif ()
  else()
    message (FATAL_ERROR "The compiler doesn't support CXX.")
  endif ()
endif ()

#==============================================================================
# CREATE nlopt_config.h
#==============================================================================

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/nlopt_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/nlopt_config.h IMMEDIATE)

# pkgconfig file
if (UNIX OR MINGW)
  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/nlopt.pc.in ${CMAKE_CURRENT_BINARY_DIR}/nlopt.pc @ONLY)
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/nlopt.pc DESTINATION ${RELATIVE_INSTALL_LIB_DIR}/pkgconfig)
endif ()

#==============================================================================
# nlopt LIBRARY TARGET (SHARED OR STATIC)
#==============================================================================

set (NLOPT_HEADERS
  src/api/nlopt.h ${PROJECT_BINARY_DIR}/nlopt.hpp ${PROJECT_BINARY_DIR}/nlopt.f
)

add_custom_command (OUTPUT nlopt.hpp
  COMMAND ${CMAKE_COMMAND} -DAPI_SOURCE_DIR=${PROJECT_SOURCE_DIR}/src/api -P ${PROJECT_SOURCE_DIR}/cmake/generate-cpp.cmake
  MAIN_DEPENDENCY src/api/nlopt-in.hpp)
add_custom_target (generate-cpp ALL DEPENDS nlopt.hpp)

add_custom_command (OUTPUT nlopt.f
  COMMAND ${CMAKE_COMMAND} -DAPI_SOURCE_DIR=${PROJECT_SOURCE_DIR}/src/api -P ${PROJECT_SOURCE_DIR}/cmake/generate-fortran.cmake
  MAIN_DEPENDENCY src/api/nlopt.h)
add_custom_target (generate-fortran ALL DEPENDS nlopt.f)

set (NLOPT_SOURCES
  src/algs/direct/DIRect.c src/algs/direct/direct_wrap.c src/algs/direct/DIRserial.c src/algs/direct/DIRsubrout.c src/algs/direct/direct-internal.h src/algs/direct/direct.�o�Do�L����S�}�H�rG����_�4��,�ظ:�6������Z4��re}�E����[�C�m��:���a��D�����e�v���9�=y��H�3�V%�o�$����?���i��5* �	�[n�����V8J݉2|M;w�Gi�'�q����V�7�S^U ��~pW�}�21MNN�;����ʼ��D�ǒeT؄�/%�y	ӯ������.5�2���� �O�&��h�� ��
{{/�ǽ�D���s��?#�����4N�i������\x��zvh� �"ɮ���9�7D3�C����U�5�I�L�����R������u�Ԧ�|�m[;I[VI�H�4��H�Ӟfz�u��ʣ��N��K�L]�)Gw�j���"��0&jp��/-�֍}	�T���#��	AJ숧�?P�kg|$�3�g-_����b%%9N�$�`B%+U+E�V���hRr���͘����@U���pP\|��WI�N8T]|{��A��2L��+Ȧ��
Cj��s�o;�Z*�K��1gE�rԱɳ���X
�-1����c���Bߪn�~���r�����۰��mN����j!��sON$[ݴ����5����=I������r������SB��W���tQ��reG]K~>k��m�t��I�A[���x��#hݐ M�/bJ5dz�BQ��cC��c�)�j��,k,�`pG�hC��U-�@J��<
c�jt	gBk����x��qb�_VӪ?��`��<|>�xv"��ѳ��#q�U�[~�lg��GB2$E�f�[���h|�ŗf|<� k:o�8&b4!�7d,�:l��hϱi�b�'�C��:t�S?�Éz��/��i\�>�Ss���<D�2�R����}�/�7E�Y"�M�Q F�XPD��* ��]N;n똽�@�V���+����h!MC�L�O�#D�'ك9���ߜ�s?����'i-L�z�F��q9')����ѡ:��郈�rKD��/�U���>�mO��Z�����6�����S�%D.w:����w��p�.ՍHq���qI���[�:�a���+�H��-dg�@L�+0	^m�>'Tb'�6yʁ���}���,M��	KL�g,��lN�<y+��d�Clޥ��$�}M\�K��w�Y��K�=A���WJ�@Rk���S1vԤJ�
�>�88Wx�9����k�}g���5Ĝ@/��0|��CSN��
AJB�,}(��7z���b�  �=>�_V���(Y'�4�DW= uB�4K
��5 ���B��i��9>��y�:Խoh��\�&�#��f?�m�#�L=� Śc�9)��}^��r�\��	/��Y"�������V��K��d��ѣ&O�V��z'��`6*++wU���ԭ��8��׵<�˵���x�Ѡ�w�%(��@�>
�u���z�����]v�QM=5�+F������q����QȮA�cb�]C�8Ω8Υ8ЩT��t1�H�z!�V�z^I+���cD��s��ևkVwE�Y/V�fu|���U������l�fgg͍��'��X\TW7c<���2�3�}��[���C����!�J|?�O'�d>���	��M���� V�"�[!��R �@�r��BY�0���bp��d�2z-�m����</��[y �� D�&�(�>�Y���˟[<dqkY)OO'Wx��`=�� ��;-��"Q�N)��������1}��{��NZ�Bj�̄��
:5�DBw���X��/�����A�w+ �bIo4���"�(�j���Xo�!���;Vf�-���͑��߹��Sp4E�x�����Ư^�߾z�ǟ�9 @p~čPa�%i@�A�����A�Rm�d2���_�'��،�Fl�%66����2 ������p;6y;�p� ��u�|�� pS�@tB��y��"���M�2�'�ɞ�>�YL3���Q��fhHv&�9կ�u�+>�_h@���xR�`6cŏ�ȥw�g��c�VJ�ܼ�b�Q�Ψ�9]�G��wP��$S��:�ԐH�o,sN�tI�� E%���&WHg�n��i�%�<v�a��`]V3o'}���������	��o��br����n4�o D�AaX����42޵Z[�LǼ�I��08�˸A>���/i�~�}9D��aܯ#�k�����9I�$�J3Ë���T���j+���kH�HWSU�א�A��SLC$:��{O2թ��+<q{�6��OL�ʧP��B�?^����s�N�»��;��/��5� ��Y���\�2WP�
��Su%��TX�V%̒��rC+�ڀ�dnK�Ҵܯί�\�ӉY��Ő�y �1Aw��J֫� f@�[;|�]n��oH?(~[qQ�������d#V)�2��������fb��	���yR�	�������A���w�#���Zn|�?L��g� �z�u]4�͆P�Ca0�0f�١j6�f#j����.��Bq9wV�3}�u�rޗx͎���L�3k���l4���F��0�Mg@b�:9��=o���#������V�����j9Չ�<�d��iMx�=̓ד�,�<�6N���?֩��[�b1�<�X��N,��$j�u?.UR�&�ݱ\Q,;�	�;�ι|n������):��5{�J{]��}��,f��K�B�YtE��+�����B�~a�#���ʆ�2���&�B�/�Gk[�CU&ni�*繘�q��i����d���	�6auM"o���Il51n(��]94S�fWKM�sXj��C+��G68��w��q����p7�#:�(���	��X�z�/��-s֕<G�yKq��B{?=eb���sS�I��e��|��̰`�!����^�w��j���W�-�%L-�U�wJ<�F�΍6z��E��ޣ+Z�2�Z�8R�/ވ���{�`�VP_���oL��mo��+�y�}3���?3��d�6�[KF�A��dX<����S{s�9�KԷ�L���w(��'��ub9ް䊣����q��{��5�����V��{�WR�O\�y��h��w��_��9|EO�<qT��ݒ'6�Zl ��XB��svL�9<)�Ɏ]О��%V��0XF���X"��/>*��iG+���J9�/�e|�T��@_��e�c�v�3MC��E�,t�E�#�v��{���|U��)����Z����wcP:>7;�/�2������$������)��ͥ��@k�m4�?3*�lHU�̂�4��x��^�ν��z�v���쓪z �,Cnj��+��*	��L�>?p+}�i�jÐv��0�3�&J1��n��2H�¼XB/�>����C���llZ^*i�Жd͚X��K��\�Y��Wt	�{����\�B��H	~e����[�}�k�0q�𵽃9��յY�Ru`d��aV��;��&.�M��R���?3^
��,[n���;�����x��x��x�w��{��{��o/}u��>�Q�|K6r�y�Ff?޽y�G&H��+�J?A*���,B�s�*:��t�I�0�{0�f��#8���$��@�,�Qz���ݮOZ�����v����U��\L�N+cM��ﳋ3ޥ�/�7������ϗ���_.�Gg�B���\��F�s�\���΋���gy���<˛�{��`
T���*�0�h���zID#`Ep9�+
���LQ����7)FbBz�����S�4빞t��Z)�@�)�+ȐB�e^z0���Z�\��m��/)/)/)��f��+y�EU
���4^�5}jF���U-RbA/�Af��H�#x���LP+Y�TB�r����%�\�������C�|�ާ��Ô��Q���`�ɖ-GP����JG��q#")�s��]�m+���Ґ��R��#޲Μ�ַ2r�T�-Yg���L�Z�"'�$�Б.)\'���.���7R
n)�t�M/"v�ڏ��'L�����&5�H�0'�O����س֓���"^9�=�ر�gƻ.9CPT_PYTHON)
  find_package (PythonInterp)
  find_package (PythonLibs ${PYTHON_VERSION_STRING} EXACT)
  find_package (NumPy)
endif ()

if (NOT DEFINED INSTALL_PYTHON_DIR AND NOT CMAKE_CROSSCOMPILING)
  execute_process (COMMAND ${PYTHON_EXECUTABLE} -c "import sysconfig, os; print(sysconfig.get_path('platlib').replace(sysconfig.get_path('data'), '').lstrip(os.path.sep))"
                    OUTPUT_VARIABLE INSTALL_PYTHON_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  file (TO_CMAKE_PATH "${INSTALL_PYTHON_DIR}" INSTALL_PYTHON_DIR)
endif ()

if (NOT DEFINED INSTALL_PYTHON_DIR)
  if (WIN32)
    set (INSTALL_PYTHON_DIR Lib/site-packages)
  else ()
    set (INSTALL_PYTHON_DIR ${INSTALL_LIB_DIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)
  endif ()
endif ()

if (NLOPT_GUILE)
  find_package (Guile)
endif ()

if (NLOPT_SWIG)
  find_package (SWIG)
  if (SWIG_FOUND)
    add_subdirectory (src/swig)
  endif ()
endif ()

if (NLOPT_OCTAVE)
  find_package (Octave)
endif ()

if (NLOPT_MATLAB)
  find_package (Matlab COMPONENTS MX_LIBRARY MAIN_PROGRAM)
endif ()

if (OCTAVE_FOUND OR Matlab_FOUND)
  add_subdirectory (src/octave)
endif ()

if (NLOPT_TESTS)
  enable_testing ()
  add_subdirectory (test)
endif ()

set (CPACK_PACKAGE_NAME          "${CMAKE_PROJECT_NAME}")
set (CPACK_PACKAGE_VERSION_MAJOR "${NLOPT_MAJOR_VERSION}")
set (CPACK_PACKAGE_VERSION_MINOR "${NLOPT_MINOR_VERSION}")
set (CPACK_PACKAGE_VERSION_PATCH "${NLOPT_BUGFIX_VERSION}")
set (CPACK_SOURCE_GENERATOR      "TBZ2;TGZ"    )
set (CPACK_BINARY_STGZ           OFF CACHE BOOL "STGZ")
set (CPACK_BINARY_TBZ2           ON  CACHE BOOL "TBZ2")
set (CPACK_BINARY_TGZ            ON  CACHE BOOL "TGZ")
set (CPACK_BINARY_TZ             OFF CACHE BOOL "TZ")
set (CPACK_SOURCE_IGNORE_FILES ".git;/build;.*~;${CPACK_SOURCE_IGNORE_FILES}")

set (CPACK_SOURCE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${NLOPT_VERSION_STRING})

include (CPack)

# configuration files
export (TARGETS ${nlopt_lib} NAMESPACE NLopt:: FILE ${PROJECT_BINARY_DIR}/NLoptLibraryDepends.cmake)

# Install the export set for use with the install-tree
install(EXPORT NLoptLibraryDepends
        NAMESPACE NLopt::
        DESTINATION ${RELATIVE_INSTALL_CMAKE_DIR}
        COMPONENT Development)

# Create a NLOPTConfig.cmake file for the use from the install tree
# and install it
set (NLOPT_LIBRARIES "NLopt::${nlopt_lib}")

set (NLOPT_CMAKE_DIR "${INSTALL_CMAKE_DIR}")
file (RELATIVE_PATH rel_include_dir "${NLOPT_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
list (APPEND RELATIVE_NLOPT_INCLUDE_DIRS ${rel_include_dir})

file (RELATIVE_PATH rel_lib_dir "${NLOPT_CMAKE_DIR}" "${INSTALL_LIB_DIR}")
list (APPEND RELATIVE_NLOPT_LIB_DIR ${rel_lib_dir})

configure_file (cmake/NLoptConfig.cmake.in NLoptConfig.cmake @ONLY)
configure_file (cmake/NLoptConfigVersion.cmake.in NLoptConfigVersion.cmake @ONLY)
install (FILES
          ${CMAKE_CURRENT_BINARY_DIR}/NLoptConfig.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/NLoptConfigVersion.cmake
         DESTINATION ${RELATIVE_INSTALL_CMAKE_DIR}
         COMPONENT Development)
