SET( evolve_SOURCE Boundary.C Element.C Evolve.C Flux.C Grid.C Initialize.C MakeHighRes.C Node.C TimeData.C)

ADD_LIBRARY(evolve ${evolve_SOURCE})

ADD_CUSTOM_TARGET( evolve_RCS
  COMMAND ci -l ${evolve_SOURCE} CMakeLists.txt
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
