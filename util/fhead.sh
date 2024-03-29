#!/bin/sh
# 
# usage is mkhead.sh ${CMAKE_CURRENT_SOURCE_DIR} where we assume that
# ${CMAKE_CURRENT_SOURCE_DIR} has subdirectories src and include and
# the function names from src/*.c are put into the include/fasp_functs.h
set +x
/bin/cat $1/src/*.c \
        | awk -v name="faspxx_functs.h" -f mkheaders.awk \
        > $1/include/faspxx_functs.h
set -x
