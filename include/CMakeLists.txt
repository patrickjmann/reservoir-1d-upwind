SET( include_SOURCE basic.h  boundary.h cloption.h  element.h  grid.h
            node.h  pexit.h  physical.h proto.h  ptime.h )
ADD_CUSTOM_TARGET( include_RCS
  COMMAND ci -l ${include_SOURCE} CMakeLists.txt
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
