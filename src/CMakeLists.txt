
include_directories (${CMAKE_SOURCE_DIR}/inc/)


add_library(inflightlib channel.cxx detector.cxx minInstance.cxx fourmomentum.cxx sterile_flux.cxx bounds.cxx decayrates.cxx)
target_include_directories (inflightlib PUBLIC ${CMAKE_SOURCE_DIR}/inc/)
target_link_libraries(inflightlib ${ROOT_LIBRARIES} cuba)



add_executable(inflight inflight.cxx)
target_link_libraries (inflight LINK_PUBLIC inflightlib)
