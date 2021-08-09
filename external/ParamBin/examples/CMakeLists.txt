
cmake_minimum_required(VERSION 3.19)

add_executable(example1 example1.cpp)
add_executable(example2 example2.cpp)
add_executable(example3 example3.cpp)
add_executable(example4 example4.cpp)
add_executable(example5 example5.cpp)
add_executable(example6 example6.cpp)

target_link_libraries(example1 PRIVATE parambin)
target_link_libraries(example2 PRIVATE parambin)
target_link_libraries(example3 PRIVATE parambin)
target_link_libraries(example4 PRIVATE parambin)
target_link_libraries(example5 PRIVATE parambin)
target_link_libraries(example6 PRIVATE parambin)

if(WIN32)
	foreach(ARG example1 example2 example3 example4 example5 example6)
		add_custom_command(TARGET ${ARG} 
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:parambin> ${CMAKE_CURRENT_BINARY_DIR}
			COMMENT "Copy parambin dll file to ${CMAKE_CURRENT_BINARY_DIR} directory" VERBATIM
		)
	endforeach()
endif()

add_custom_command(TARGET example1 POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy
					${CMAKE_CURRENT_SOURCE_DIR}/input1.txt
					${CMAKE_CURRENT_BINARY_DIR}/input1.txt
)

add_custom_command(TARGET example1 POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy
					${CMAKE_CURRENT_SOURCE_DIR}/input2.txt
					${CMAKE_CURRENT_BINARY_DIR}/input2.txt
)

add_custom_command(TARGET example1 POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy
					${CMAKE_CURRENT_SOURCE_DIR}/argon.txt
					${CMAKE_CURRENT_BINARY_DIR}/argon.txt
)

