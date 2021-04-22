configure_file("modulefiles/anamerger" anamerger @ONLY)
configure_file("setup.sh.in" setup.sh @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/anamerger DESTINATION share/modulefiles/)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.sh DESTINATION share/)
