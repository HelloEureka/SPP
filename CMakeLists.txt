cmake_minimum_required(VERSION 3.17)
project(SPP)

# 添加math库
LINK_LIBRARIES(m)

#set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_STANDARD 11)

add_executable(SPP main.c com.c SatPos.c CMat.c rnx.c)
