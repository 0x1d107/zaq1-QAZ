cmake_minimum_required(VERSION 3.22)
project(inv_q_filtr)

set(ENABLE_PRECOMPILED_HEADERS OFF)

add_subdirectory(source)
add_subdirectory(scripts)
add_compile_options(-fsanitize=address)
add_link_options(-fsanitize=address)
