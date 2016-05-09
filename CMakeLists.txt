cmake_minimum_required (VERSION 2.8)

set(PROJECT VectorFEM)

project (${PROJECT})

# -Wall     - отображение всех предупреждений
# -O0       - отключает оптимизацию. “олько переменные, объ€вленные register, сохран€ютс€ в регистрах
# -pedantic - предупреждени€ компил€тора (св€зано с -Wall)
# 
# add_definitions - указание флагов компил€ции дл€ текущей директории и поддиректорий
if (MSVC)
	add_definitions (/W4) # максимальный уровень предупреждений /W4 дл€ компил€тора от Microsoft
elseif (CMAKE_COMPILER_IS_GNUCXX)
	add_definitions (-Wall -pedantic)
else ()
	message ("Unknown compiler")
endif ()

set (SOURCE_ROOT src)

aux_source_directory( ${SOURCE_ROOT} SOURCES )
file(GLOB HEADERS "${SOURCE_ROOT}/*.h")

include_directories(src)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

add_executable (${PROJECT} ${HEADERS} ${SOURCES})
