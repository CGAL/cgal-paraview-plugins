set(BUILD_SHARED_LIBS TRUE)

#Paraview's macro to add the plugin. It takes care of all the vtk
#and paraview parts of the process, like link and integration
#in the UI
if(ParaView_FOUND AND PARAVIEW_USE_FILE)
  ## ParaView 5.6
    add_subdirectory(IsotropicRemeshingFilters)
else()
  ## ParaView 6.x
  vtk_module_scan(
    MODULE_FILES      "${CMAKE_CURRENT_SOURCE_DIR}/IsotropicRemeshingFilters/vtk.module"
    REQUEST_MODULES   IsotropicRemeshingFilters
    PROVIDES_MODULES  provided_modules)
  
  vtk_module_build(MODULES ${provided_modules})
  
  paraview_add_plugin(IsotropicRemeshingFilter
    VERSION "1.0"
    MODULES IsotropicRemeshingFilters)
endif()

