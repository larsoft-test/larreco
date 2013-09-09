
cet_add_compiler_flags( CXX -Werror )

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
	                ${ROOT_CORE}
                        ${ROOT_CINT} 
                        ${ROOT_RIO}
                        ${ROOT_NET}
                        ${ROOT_HIST} 
                        ${ROOT_GRAF}
                        ${ROOT_GRAF3D}
                        ${ROOT_GPAD}
                        ${ROOT_TREE}
                        ${ROOT_RINT}
                        ${ROOT_POSTSCRIPT}
                        ${ROOT_MATRIX}
                        ${ROOT_PHYSICS}
                        ${ROOT_MATHCORE}
                        ${ROOT_THREAD}
                        ${ROOT_MINUIT}
        )
#                        ${ROOT_GEOM}
#                        ${ROOT_XMLIO}
#                        ${ROOT_GDML}

install_headers()
install_fhicl()
install_source()
