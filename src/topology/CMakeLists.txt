set(main_target "topology")

add_executable(${main_target} main.cpp)
target_sources(${main_target}
  PRIVATE
    main.cpp
)

target_compile_options(${main_target} PUBLIC ${compile_options_tmp})

target_link_libraries(${main_target} PRIVATE
  STOKES
  ${TP_LIBRARY}
  ${MKL_LIBRARIES}
  ${NLopt_LIBRARY}
  iomp5
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})