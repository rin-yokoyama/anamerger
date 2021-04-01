configure_file( config_gen.sh.in ${CMAKE_CURRENT_BINARY_DIR}/generated/config_gen.sh @ONLY)

install(FILES anamerger_config.yaml DESTINATION config)
install(FILES configallnuclides.sh DESTINATION config)
install(FILES config_gen.py DESTINATION config)  
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/generated/config_gen.sh DESTINATION config)
