# Find Cell Broadband Enginge (CBE)
# Varaibles
#  CBE_FOUND
#  CBE_INCLUDE_DIR
#  CBE_LIBRARY
#  SPU_CC
#  SPU_CXX
#  PPU_CC
#  PPU_CXX
#  PPU_EMBEDSPU
# Macro
#  CBE_COMPILE_C (
#      GENERATED_CBE_FILES[out] 
#      SPU_C_SOURCES[in] 
#      SPU_C_FLAGS[in]
#      PPU_EMBED_FLAGS[in]
#      PPU_C_SOURCES[in] 
#      PPU_C_FLAGS[in]
#  )

find_path( CBE_SPE_INCLUDE_DIR libspe2.h )

find_library ( CBE_SPE_LIBRARY spe2 )

set(CBE_INCLUDE_DIR ${CBE_SPE_INCLUDE_DIR})
set(CBE_LIBRARY ${CBE_SPE_LIBRARY})

if ( NOT SPU_CC )
    find_program( SPU_CC NAMES spu-gcc43 )
endif ( NOT SPU_CC )

if ( NOT SPU_CXX )
    find_program( SPU_CXX NAMES spu-g++43 )
endif ( NOT SPU_CXX )

if ( NOT PPU_CC )
    find_program( PPU_CC NAMES ppu-gcc43 )
endif ( NOT PPU_CC )

if ( NOT PPU_CXX )
    find_program( PPU_CXX NAMES ppu-g++43 )
endif ( NOT PPU_CXX )

#if ( NOT PPU_AR )
#    find_program( PPU_AR NAMES ppu-ar )
#endif ( NOT PPU_AR )

if ( NOT PPU_LD )
    find_program( PPU_LD NAMES ppu-ld )
endif ( NOT PPU_LD )


if ( NOT PPU_EMBEDSPU )
    find_program( PPU_EMBEDSPU NAMES ppu-embedspu )
endif ( NOT PPU_EMBEDSPU )

# handle the QUIETLY and REQUIRED arguments and set CBE_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CBE DEFAULT_MSG
  PPU_CC
  SPU_CC
  PPU_CXX
  SPU_CXX
  PPU_EMBEDSPU
  #PPU_AR
  #SPU_AR
  #PPU_LD
  #SPU_LD
  #own
  CBE_LIBRARY
  CBE_INCLUDE_DIR
  )
  
MARK_AS_ADVANCED(PPU_CXX SPU_CXX PPU_CC SPU_CC PPU_EMBEDSPU PPU_LD)

MARK_AS_ADVANCED(CBE_LIBRARY 
                 CBE_INCLUDE_DIR 
                 CBE_SPE_INCLUDE_DIR 
                 CBE_SPE_LIBRARY)

#generates the spu side. returns the ppu-embedspu files. 
#expects all spu files in <current_dir>/CellKernel/<spu-source-files>
macro(COMPILE_SPU GENERATED_SPU_FILES)
    foreach( SPU_SOURCE_NAME ${SPU_C_SOURCES} )  
        set (SPU_OBJ ${SPU_SOURCE_NAME} )
        set (SPU_SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/CellKernel/${SPU_SOURCE_NAME}.c )
        
        add_custom_command(
            OUTPUT ${SPU_OBJ}.o
            COMMAND ${SPU_CC} 
            ARGS ${SPU_C_FLAGS} -o ${SPU_OBJ}.o ${SPU_SOURCE_FILE}
            DEPENDS ${SPU_SOURCE_FILE}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        ) 
        
        add_custom_command(
            OUTPUT ${SPU_SOURCE_NAME}-embed.o
            COMMAND ${PPU_EMBEDSPU} 
            ARGS ${PPU_EMBED_FLAGS} ${SPU_SOURCE_NAME} ${SPU_SOURCE_NAME}.o ${SPU_SOURCE_NAME}-embed.o
            DEPENDS ${SPU_OBJ}.o
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        ) 
        
        set( SPU_EMBEDS ${SPU_EMBEDS} ${SPU_OBJ}-embed.o)
    endforeach()
    
    set (GENERATED_SPU_FILES ${SPU_EMBEDS})
endmacro(COMPILE_SPU)