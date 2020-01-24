project(target)

add_subdirectory(solver)
add_subdirectory(utils)

add_executable(main main.cpp)

target_sources(
    main
    PRIVATE
        main.cpp
    PUBLIC
)

target_link_libraries(
    main
    solverlib
    utils
)
