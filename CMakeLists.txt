cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(CGAL_Isotropic_remeshing_filter)
#If the plugin is used internally (inside Paraview's source directory),
#then we don't need to call find_package.
if (NOT ParaView_BINARY_DIR)
  find_package(ParaView REQUIRED)
endif()
#Find CGAL
find_package(CGAL)
if(CGAL_FOUND)
  include( ${CGAL_USE_FILE} )
endif(CGAL_FOUND)
if(ParaView_FOUND)
  #If the plugin is used internally we don't need to include.
  if (NOT ParaView_BINARY_DIR)
      include(${PARAVIEW_USE_FILE})
  endif(NOT ParaView_BINARY_DIR)
  
  #Paraview's macro to add the plugin. It takes care of all the vtk 
  #and paraview parts of the process, like link and integration
  #in the UI
  ADD_PARAVIEW_PLUGIN(IsotropicRemeshingFilter "1.0"
    SERVER_MANAGER_XML IsotropicRemeshingFilter.xml
    SERVER_MANAGER_SOURCES vtkIsotropicRemeshingFilter.cxx)
else()
  message("WARNING : The Paraview plugins need Paraview, so they won't be compiled.")
endif(ParaView_FOUND)
# Link with CGAL
if(CGAL_FOUND)
  if(ParaView_FOUND)
    target_link_libraries(IsotropicRemeshingFilter LINK_PRIVATE
      CGAL::CGAL 
      ${Boost_LIBRARIES})
  endif(ParaView_FOUND)
endif (CGAL_FOUND)
