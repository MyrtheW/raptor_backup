cmake_minimum_required (VERSION 3.18)

add_library ("raptor_upgrade" STATIC
        raptor_upgrade.cpp
        load_hibf.cpp
        get_fpr.cpp
        )
target_link_libraries ("raptor_upgrade" PUBLIC "raptor_interface")
