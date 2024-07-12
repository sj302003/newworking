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
  src/algs/direct/DIRect.c src/algs/direct/direct_wrap.c src/algs/direct/DIRserial.c src/algs/direct/DIRsubrout.c src/algs/direct/direct-internal.h src/algs/direct/direct.÷oæDo‚LË¯ûÏS∏}§HærGâÏ·â¡_°4¨€,˘ÿ∏:ü6¥øîíÒ Z4≤ære}µE÷⁄Òú[¥C∏mΩø:üè˘a –DΩ˙¨´©e¶vˇ‚Í9»=yêÔôHÅ3îV%Ÿo£$¥Ö≤‡?¡Öi¬ììÒß5*  	…[näîë–ÜV8J›â2|M;wûGiñ'Î¢q•®ê˛VΩ7›S^U ô‡~pW±}Í21MNNæ;¯Ó‰‰ ºÆ‡D·«íeTÿÑ±/%åy	”Ø∑ÎÈõûÇ.5¢2Î¸˝Ú ˇOŒ&ãh—ø Ωô
{{/É«ΩåD◊ÀÏâsëø?#‡îü©˜4NÍiåØë∫„¯\x˝Ízvh⁄ œ"…Æã˜ç9¿7D3Ûõ∏C´Ñˇ‚U™5¡IüLß˜À˚ÃRÉ±˘≥ù¢u€‘¶ç|÷m[;I[VIåHÂ4–ŸHø”ûfz∞u‡à £ÅÑN®ÀK¢L]ù)GwäjΩ…˚"π∆0&jpÎÊ/-’÷ç}	‘T™ıÛ#—Ø	AJÏàßÒ?PÖkg|$«3˛g-_¢¿ì˚b%%9N¥$ˇ`B%+U+EîVÇàÄhRr¶´ªÕòü®‡Ê@Uü•öpP\|çÖWI÷N8T]|{âìAÛÚ2LÀå+»¶ëê
Cj®êsÂo;≥Z*òKÖê1gEÖr‘±…≥çø€X
∂-1¡£§ècÍ¸‰Bﬂ™nñ~¿É˙r”¥¸∫Â€∞©ˆmN®¿˝•j!ÅÀsON$[›¥∂ù˙∆5ÔÎ¯∆=IÉ„µ¡™ørÊ‡å¬‰†ÿSBø«WÙ˜átQ·„≥reG]K~>kõèmÎt∏çIöA[åøx®≤#h›ê Mª/bJ5dz¥BQÂ”cCö“c›)Öj™¬,k,Ä`pGáhCµU-±@J≥±<
cêjt	gBkìÂÛ xÍ°qb¡_V”™?íÚ`∫ö<|>…xv"π”—≥º∞#qôUÈ[~„©lgÜÇGB2$Eàf°[∫ëıh|Ã≈óf|<å k:oÒ8&b4!ç7d,ç:lå†hœ±iûbá'∫Cˆ’:tπS?ú√âz≥Ä/ÍØêi\Æ>Ss´èπ<DÓ2„R™©ÁÙ}Ω/§7E‘Y"¯MœQ FÄXPDÄ¢* õº]N;nÎòΩÑ@ÆVÍ†„+ûÒªüÊÈh!MCÏªLÛOú#D¸'ŸÉ9ˆÁ≤Ëﬂú»s? ÛÇ≥'i-LèzáFú¬q9')¯∞˘Òå—°:ÑêÈÉàÑrKD Ó/ÎUëà◊>ÊmO´ÉZ⁄ÈîÏÍ‚6Ÿ»¸ó¡Sæ%D.w:Ê§Ò›≈w¯⁄pí.’çHqãé≈qIäö´[í:æaµöÛ+‹Hø -dg¥@L»+0	^m‹>'Tb'¯6y Å∏ÅÒí}≥îª,MåÍé	KL÷g,ølNã<y+ÈÊdÛClﬁ•ºÇ$Ñ}M\ãKÅÆwÖY™˝Kò=AˇÛÑÔWJÎ@Rk™≥ïS1v‘§Jı
π>˛88Wxó9≥µÃ€kÂ}gˆÈ”5ƒú@/ÃÔ0|â¡CSNÊ⁄
AJBô,}(±≠7z¥•´b“  √=>ËÖ_Vè®±(Y'Ÿ4©DW= uB∂4K
˛§5 πÃ≈Bá◊i·Ï9>ü˙yÈ:‘Ωoh˝’\‡&Á#Öüf?≠mÏ#íL=˘ ≈öcﬂ9)˚”}^ò◊rª\ﬁ‡	/¶´Y"È˝ç∂«Â¡VÍÓK≥Ádé”—£&O·V£Æz'Öı`6*++wU⁄œÙ‘≠ø‚ª8¡˙◊µ<‰ÀµÀÒé˚x∑—†Éwî%(¬˜@˚>
÷uëÚÈzÎµ „áû]v…QM=5÷+FÚêÚ¯Å¡ qàÜ¸ÊQ»ÆA»cbê]Cü8Œ©8Œ•8–©Täèt1Ë¶Hùz!òV“z^I+∂•ácD⁄às¶À÷ákVwE∫Y/V¡fu|ú‡°ÓU†∂˜˛ØÖlöfggÕç‚∏ö'≈ˆX\TW7c<óå˜2ñ3æ}£Óº[’¯⁄C¶À‡â!¨J|?öO'ãd>ô™Õ	”“Mà·÷“ V›"ª[!’‰R ∫@ûrâ¶BYÿ0˚Ú›bpü·õdò2z-°mÉ˙Æà</ò©[y £˜ D˜&Ä(≠>ÔüY˘Í‘Àü[<dqkY)OO'Wx©æ`=∂≈ ÛÏ;-ı≥"QﬁN)Ñ∞üâÛÊœÙ1}ª˙{Ò“NZÑBjüÃÑ„÷
:5¡DBwµ‡◊X’/ˇº´º•A™w+ ˙bIo4¥†µ"˘( j˚›ÚXoπ!æì∆;Vfù-¿«ÕëÎ‚ﬂπ‚≈Sp4E“x˝ÍÕœÔ∆Ø^åﬂæz˚«ü±9 @p~ƒçPaû%i@ÈA∞‹ÃÔ≥ÈA≠Rm“d2õçì_Ô'ã¥ÿåóFl≤%66ÖÄ¬Ï2 Á˙¿É±Ñp;6y;¡pÙ Ö„u≤|Û› pS≈@tB˝·yïË"ª¥ÌM‹2’'π…ûÈÉ>˝YL3“˙§Q∂˚fhHv&≈9’ØáuÚ+>„_h@ØŸ‰xR `6c≈è‡»•wÂg≤¸c≤VJ§‹ºêbÉQ¿Œ®æ9]∆GÉ»wPì€$S·∂:î‘êHÑo,sN∫tIı‰ç E%ô•‡&WHg∞n∫ÜiË%Ñ<vÊa©ƒ`]V3o'}¶ôú∑†˘ƒ˛…	˛ÇoÇ˘br≠ñ˝æn4Ào DôAaXø⁄⁄÷42ﬁµZ[ÚêL«ºÅI‡Ñ08ïÀ∏A>…¿·/iÒ~¸}9DÎÙa‹Ø#∂kö¡Ñ†⁄9Iﬂ$§J3√ãÑ´ÍåTπãçj+≥≤äkH˛HWSUÅ◊ê»Aí‹SLC$:ﬁÎ{O2’©‘›+<q{Ä6‚œOLƒ ßPòÕB¸?^ÓÃ¬ÊsÌN˘¬ª‹≈;úà/ΩÛ5© —‚Y≈≥–\”2WPÄ
ßòSu%Ó»TX–V%ÃíïßrC+Ò⁄ÄÆdnKß“¥‹ØŒØÎ\∆”âY‰≈≈ê©y Ø1AwÉˇJ÷´ä f@ò[;|π]n˘ÊoH?(~[qQ±äèπóèπd#V)ì2úâøÔ»ƒ‹fbæç	˙≠ÂyR´	ÉŸ¸öÉﬁŒAÅówâ#¯„ûÁZn| ?L˛ÂgŸ æz˚u]4ÓÄÕÜPëCa0∞0f∏Ÿ°j6≤f#jêæÅß.ËŒBq9wVﬁ3}íuÂ≠rﬁóxÕéøºıL◊3kÁ¡Äl4äç≈F£ÿ0≥Mg@bŒ:9Å∫=oô·À#áˇ†É“ÁVä˘í¡Èj9’âè<õdçÏ°iMxÄ=ÕÉ◊ìÃ,‰<ç6N£ç∑?÷©Öë[∏b1£<âX—»N,ØÉ$j•u?.UR¨&æ›±\Q,;‰	 ;á·ææ|n»ØâßÍ‰):ùæ5{ÔJ{]ë¢}∞ó,f»ÍK≤BˆYtEè†+¶ÚÚ“∫B¶~aÒ#©¿¸ Üí2˙í≠&’BÓ/©Gk[ÈCU&niË*ÁπòÂq˘¿ÔúäiÀÀùÂdì∞Æ	Û6auM"oì®ÆIl51n(áõ]94SıfWKM∂sXj≤ùC+Ü∞G68ıÂw◊˙q†Úª˚Ûp7ç#:â(˛üˇ	¨íXÊzÁ/ÔÙ-s÷ï<GÊyKqπ≈B{?=eb⁄˛¿sSÍIµîe¸»|°òÃ∞`°!É∫Í»^ˇw∫¬jü™ÿW§-å%L-ÛU„wJ<òF—Œç6zç˛E“Ëéﬁ£+ZÍ2∆Z¢8R≈/ﬁà‚ »{˝`‚ÉVP_ﬂˆ’oL¯·moÍ‡´+ÜyÃ}3™Éˇ?3¸ØdÃ6Ü[KFÏA√Ú∏dX<ó‚ËÍS{s…9≥K‘∑’Lâ˙æw(‹Ô'·ÿub9ﬁ∞‰ä£ﬂÃáèq≈¸{£€5£´¿ö’V’{˘WRΩO\¯y«Ùh≈¯w˜ƒ_¿øí9|EOñ<qTÚƒ›í'6ﬁZl ÚXBÇµsvL 9<)Í…é]–ûÏûı%V‡„0XF¥ú£X"ä•/>*∏åiG+¿ºﬂJ9Ò/¸e|‹Të¯@_∂—eÍcõvë3MCıÒEß,tÒEÚ#èvëÇ{êº¯|Uü˙)Ì˚®¶Zè˙æ≠wcP:>7;¶/›2Ì›Ï›¯ó$€¿ñâÔ˚)ô˙Õ•«˘@kÒm4Ë?3*√lHUêÃÇ‹4âÒxê˙^ŸŒΩÜ€zçvÈı°Ïì™z É,Cnj⁄¶+®Ã*	¿ÌL¬>?p+}ıiºj√êvŸ 0∫3î&J1∏ón•˜2H≈¬ºXB/ñ>≈„å∆C´‹òllZ^*iß–ñdÕöXˇªK¯Ì\€Y°ŸWt	€{∑ı˙’\˚BóH	~eóÄπÚæä[∞}ëk∞0q˜µΩÉ9ı»’µY¡Ru`dïπaVËƒ;™∞&.≤M†ÖR¥¸Ÿ?3^
Â∑,[n¸¯«;«œÏÁ˜xÈ˜xÈ˜xÈwó{ºÙ{ºÙo/}u˜†>ùQõ|K6r˛yÏFf?ﬁΩyÖG&H·ÿ+îJ?A*îÎ˘,B‚sÛ*:Ê’t∏IÔ0•{0…f¶€#8Ñö…$ÔÍ@ú,˝Qzª¥æ›ÆOZ üœÒ‰v≤Ä⁄‡UöΩ\LÚ¢°N+cMûõÔ≥ã3ﬁ•Ô≥/≈7È˜Ÿ›Ù˝œóµŸ˝_.ΩGg…Bü«Ù\∞ÄFòsü\˜˙ëŒã¥∏ÁgyÅ≥â<Àõ„{êù`
TÔÇÊÃ*±0èhª¨ˇzID#`Ep9È+
†¸©LQ˜≤¶Œ7)FbBzáàÌá°S¶4Îπût†ÍπZ)@ÿ)¿+»êBäe^z0µãÃZç\ÄÜm⁄¿/)/)/)ö∞fØ÷+y“EU
™ë°4^ª5}jF¿“‘U-RbA/ÕAfı˚H†#xÌ⁄˘LP+Y˜TBüràìÆ∏%ø\Øë¯Û√µﬂC∆|Ôﬁß÷Û√îüŒQèÒòÚ`Ûè…ñ-GP∑†˜µJG¿Õq#")˘sÊµ˚]Œm+º¡ñ“êÛ›R®ﬁ#ﬁ≤ŒúÚ÷∑2r¸T¡-YgÚÎﬁLåZ∂"'›$Îµ–ë.)\'øûœ.Üó≤7R
n)ætË∫M/"vÈ⁄èÄŸ'L™é«ÍÍ&5πH’0'≈Oßº¿Ëÿ≥÷ì¥∂º"^9ª=âÿ±æg∆ª.9CPT_PYTHON)
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
