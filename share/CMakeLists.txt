configure_file("modulefiles/anamerger" anamerger @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/anamerger DESTINATION share/modulefiles/)