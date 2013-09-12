
art_make( 
          LIB_LIBRARIES HarrisVertexFinder_module
                        ${BACKTRACKER_SERVICE}
                        ${RECOBASE}
                        ${SIMULATION}
                        ${FILTERS}
                        ${GEOMETRY}
                        ${GEOMETRY_SERVICE}
                        ${RECOOBJECTS}
                        ${ANALYSISBASE}
                        ${LARPROPERTIES_SERVICE}
                        ${DETECTORPROPERTIES_SERVICE}
                        ${UTILITIES}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${CLHEP}
			${ROOT_BASIC_LIB_LIST}
                        ${ROOT_MINUIT}
        )

install_headers()
install_fhicl()
install_source()
