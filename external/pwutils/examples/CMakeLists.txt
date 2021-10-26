
cmake_minimum_required(VERSION 3.19)

add_executable(report_example report_example.cpp)
add_executable(read_example read_example.cpp)
add_executable(indexing_example indexing_example.cpp)

target_link_libraries(report_example PUBLIC pwutils)
target_link_libraries(read_example PUBLIC pwutils)
target_link_libraries(indexing_example PUBLIC pwutils)

if(WIN32)
    foreach(ARG report_example read_example indexing_example)
        add_custom_command(TARGET ${ARG}
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:pwutils> ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Copy pwutils dll file to ${CMAKE_CURRENT_BINARY_DIR} directory" VERBATIM
        )
        add_custom_command(TARGET ${ARG}
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:json11> ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Copy json11 dll file to ${CMAKE_CURRENT_BINARY_DIR} directory" VERBATIM
        )
    endforeach()
endif()

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})



