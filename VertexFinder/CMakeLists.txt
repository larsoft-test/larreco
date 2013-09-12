
cet_remove_compiler_flag( CXX -Werror )

art_make( BASENAME_ONLY
          EXCLUDE FeatureVertexFinder_module.cc
          MODULE_LIBRARIES 
                        ${BACKTRACKER_SERVICE}
                        ${FILTERS}
                        ${RECOBASE}
                        ${GEOMETRY}
                        ${GEOMETRY_SERVICE}
                        ${LARPROPERTIES_SERVICE}
                        ${DETECTORPROPERTIES_SERVICE}
                        ${UTILITIES}
                        ${SIMULATION}
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
         )


simple_plugin ( FeatureVertexFinder "module"
                        RecoAlg
                        ${BACKTRACKER_SERVICE}
                        ${FILTERS}
                        ${RECOBASE}
                        ${GEOMETRY}
                        ${GEOMETRY_SERVICE}
                        ${LARPROPERTIES_SERVICE}
                        ${DETECTORPROPERTIES_SERVICE}
                        ${UTILITIES}
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
                BASENAME_ONLY 
              )

install_headers()
install_fhicl()
install_source()


