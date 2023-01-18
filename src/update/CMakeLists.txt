cmake_minimum_required (VERSION 3.15)

add_library ("raptor_update" STATIC
        raptor_update.cpp
        load_hibf.cpp
        insertions.cpp
        rebuild.cpp
        )
target_link_libraries ("raptor_update" PUBLIC "raptor_interface")

target_link_libraries ("raptor_update" PUBLIC "raptor_interface" "chopper_layout_lib" xxHash::xxhash)