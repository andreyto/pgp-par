find_package( EXPAT REQUIRED )

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${EXPAT_INCLUDE_DIRS})

set(SOURCES base64.c BN.c BuildMS2DB.c ChargeState.c CMemLeak.c Errors.c ExonGraphAlign.c 
	FreeMod.c IonScoring.c 
	LDA.c main.c Mods.c MS2DB.c ParentMass.c ParseInput.c ParseXML.c PValue.c 
	Run.c Score.c Scorpion.c SNP.c Spectrum.c Spliced.c 
	SpliceDB.c SpliceScan.c SVM.c Tagger.c Trie.c Utils.c)

add_executable(inspect ${SOURCES})
target_link_libraries(inspect ${EXPAT_LIBRARIES} m)

install(TARGETS inspect RUNTIME DESTINATION ${PGP_INSTALL_PREFIX})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ DESTINATION ${PGP_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)


