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

  # Link with CGAL
  target_link_libraries(IsotropicRemeshingFilter PRIVATE
    CGAL::CGAL
    ${Boost_LIBRARIES})

  target_compile_definitions(IsotropicRemeshingFilter PRIVATE
    ISOTROPICREMESHINGFILTERS_EXPORT=
    CGAL_ISOTROPIC_REMESHING_NO_MODULE)
else()
set(classes
  vtkIsotropicRemeshingFilter)

vtk_module_add_module(IsotropicRemeshingFilters
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  IsotropicRemeshingFilter.xml)
endif()
