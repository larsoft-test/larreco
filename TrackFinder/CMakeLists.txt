
#cet_add_compiler_flags( CXX -Werror )

include_directories ( $ENV{GEANT4_FQ_DIR}/include )
include_directories ( $ENV{XERCES_C_INC} )

set( G4_LIB_LIST ${XERCESC}
	                ${G4FR}
	                ${G4GMOCREN}
	                ${G4OPENGL}
	                ${G4RAYTRACER}
	                ${G4TREE}
	                ${G4VRML}
	                ${G4ANALYSIS}
	                ${G4DIGITS_HITS}
	                ${G4ERROR_PROPAGATION}
	                ${G4EVENT}
	                ${G4GEOMETRY}
	                ${G4GL2PS}
	                ${G4GLOBAL}
	                ${G4GRAPHICS_REPS}
	                ${G4INTERCOMS}
	                ${G4INTERFACES}
	                ${G4MATERIALS}
	                ${G4MODELING}
	                ${G4PARMODELS}
	                ${G4PARTICLES}
	                ${G4PERSISTENCY}
	                ${G4PHYSICSLISTS}
	                ${G4PROCESSES}
	                ${G4READOUT}
	                ${G4RUN}
	                ${G4TRACK}
	                ${G4TRACKING}
	                ${G4VISHEPREP}
	                ${G4VISXXX}
	                ${G4VIS_MANAGEMENT}
	                ${G4ZLIB}
   )

art_make( BASENAME_ONLY 
          MODULE_LIBRARIES Genfit RecoAlg
                        ${SIMULATION}
                           ${RECOOBJECTS}
                           ${RECOBASE}
	                ${GEOMETRY}
	                ${GEOMETRY_SERVICE}
                        ${ANALYSISALG}
                        ${ANALYSISBASE}
                        ${UTILITIES}
                        ${MAGNETICFIELD_SERVICE}
                           ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                           ${CLHEP}
                        ${ROOT_LIB_LIST}
                        ${G4_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()


