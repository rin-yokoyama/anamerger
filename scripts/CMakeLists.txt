configure_file( run_anamerger.pbs.in ${CMAKE_CURRENT_BINARY_DIR}/generated/run_anamerger.pbs @ONLY)

install(FILES run_anamerger.py DESTINATION scripts)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/generated/run_anamerger.pbs DESTINATION scripts)
