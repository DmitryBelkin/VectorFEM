cmake_minimum_required (VERSION 2.8)

set(PROJECT VectorFEM)

project (${PROJECT})

# -Wall     - ����������� ���� ��������������
# -O0       - ��������� �����������. ������ ����������, ����������� register, ����������� � ���������
# -pedantic - �������������� ����������� (������� � -Wall)
# 
# add_definitions - �������� ������ ���������� ��� ������� ���������� � �������������
if (MSVC)
	add_definitions (/W4) # ������������ ������� �������������� /W4 ��� ����������� �� Microsoft
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
