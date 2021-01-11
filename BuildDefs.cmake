set(DATA_TYPES_MSG "Template Types Set:")

if (BUILD_VIRTUAL)
	list(APPEND DATA_DEFS "BUILD_VIRTUAL")
	set(DATA_TYPES_MSG "${DATA_TYPES_MSG} virtual;")
endif ()

# Check what Templates to initialize
string(TOLOWER ${DATA_TYPES} DATA_TYPES)
separate_arguments(DATA_TYPES NATIVE_COMMAND "${DATA_TYPES}")
string(REPLACE "," ";" DATA_TYPES "${DATA_TYPES}")
set(DATA_TYPES_LIST "${DATA_TYPES}" CACHE INTERNAL "List of template types to be built")
if ("all" IN_LIST DATA_TYPES_LIST OR
		"single" IN_LIST DATA_TYPES_LIST OR
		"float" IN_LIST DATA_TYPES_LIST)
	set(DATA_TYPES_MSG "${DATA_TYPES_MSG} float;")
	list(APPEND DATA_DEFS "BUILD_FLOAT")
	if ("all" IN_LIST DATA_TYPES_LIST OR
			"complex" IN_LIST DATA_TYPES_LIST)
		set(DATA_TYPES_MSG "${DATA_TYPES_MSG} complex<float>;")
		list(APPEND DATA_DEFS "BUILD_CFLOAT")
	endif ()
elseif ("complexfloat" IN_LIST DATA_TYPES_LIST OR
		"cfloat" IN_LIST DATA_TYPES_LIST)
	set(DATA_TYPES_MSG "${DATA_TYPES_MSG} complex<float>;")
	list(APPEND DATA_DEFS "BUILD_CFLOAT")
endif ()
if ("all" IN_LIST DATA_TYPES_LIST OR
		"double" IN_LIST DATA_TYPES_LIST)
	set(DATA_TYPES_MSG "${DATA_TYPES_MSG} double;")
	list(APPEND DATA_DEFS "BUILD_DOUBLE")
	if ("all" IN_LIST DATA_TYPES_LIST OR
			"complex" IN_LIST DATA_TYPES_LIST)
		set(DATA_TYPES_MSG "${DATA_TYPES_MSG} complex<double>;")
		list(APPEND DATA_DEFS "BUILD_CDOUBLE")
	endif ()
elseif ("complexdouble" IN_LIST DATA_TYPES_LIST OR
		"cdouble" IN_LIST DATA_TYPES_LIST)
	set(DATA_TYPES_MSG "${DATA_TYPES_MSG} complex<double>;")
	list(APPEND DATA_DEFS "BUILD_CDOUBLE")
endif ()
message(${DATA_TYPES_MSG})
set(DATA_DEFS_LIST "${DATA_DEFS}" CACHE INTERNAL "List of template definition flags")