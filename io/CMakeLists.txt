# Generally best to set the source files and use them.
# NOTE: this ends up as a semi-colon delimited list, but ADD_LIBRARY(..)
#   seems to understand it.

SET( io_SOURCE cloption.C Element.C ErrorExit.C GetId.C Grid.C MakeOutFiles.C Node.C
     NormalExit.C )

ADD_LIBRARY(io ${io_SOURCE} )

ADD_CUSTOM_TARGET( io_RCS
  COMMAND ci -l ${io_SOURCE} CMakeLists.txt
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
