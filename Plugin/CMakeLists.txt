
set(BUILD_SHARED_LIBS TRUE)

#Find CGAL
find_package(CGAL REQUIRED)
include( ${CGAL_USE_FILE} )

#Paraview's macro to add the plugin. It takes care of all the vtk
#and paraview parts of the process, like link and integration
#in the UI
if(ParaView_FOUND AND PARAVIEW_USE_FILE)
  ## ParaView 5.6
  add_paraview_plugin(IsotropicRemeshingFilter "1.0"
    SERVER_MANAGER_XML IsotropicRemeshingFilter.xml
    SERVER_MANAGER_SOURCES vtkIsotropicRemeshingFilter.cxx)
else()
  ## ParaView 6.x
  paraview_add_plugin(IsotropicRemeshingFilter
    VERSION "1.0"
    REQUIRED_ON_SERVER
    SERVER_MANAGER_XML IsotropicRemeshingFilter.xml
    SOURCES vtkIsotropicRemeshingFilter.cxx)
endif()

# Link with CGAL
target_link_libraries(IsotropicRemeshingFilter PRIVATE
  CGAL::CGAL
  ${Boost_LIBRARIES})
