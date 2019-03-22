#include <ostream>
#include <sstream>
#include "vtkInformationVector.h"
#include "vtkIsotropicRemeshingFilter.h"
#include "vtkPolyData.h"
#include "vtkInformation.h"
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

// Declare the plugin
vtkStandardNewMacro(vtkIsotropicRemeshingFilter);

namespace PMP = CGAL::Polygon_mesh_processing;
// Useful typedefs
typedef CGAL::Simple_cartesian<double>                            K;
typedef CGAL::Surface_mesh<K::Point_3>                            SM;
typedef boost::property_map<SM, CGAL::vertex_point_t>::type       VPMap;
typedef boost::property_map_value<SM, CGAL::vertex_point_t>::type Point_3;
typedef boost::graph_traits<SM>::vertex_descriptor                vertex_descriptor;
typedef boost::graph_traits<SM>::edge_descriptor                  edge_descriptor;
typedef boost::graph_traits<SM>::face_descriptor                  face_descriptor;
typedef boost::graph_traits<SM>::halfedge_descriptor              halfedge_descriptor;

// -----------------------------------------------------------------------------
// Constructor
// Fills the number of input and output objects.
// Initializes the members that need it.
vtkIsotropicRemeshingFilter::vtkIsotropicRemeshingFilter()
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

// ----------------------------------------------------------------------------
// Gets the input
// Creates CGAL::Surface_mesh from vtkPolydata
// Calls the CGAL::isotropic_remeshing algorithm
// Fills the output vtkPolyData from the result.
int vtkIsotropicRemeshingFilter::RequestData(vtkInformation *,
                                             vtkInformationVector **inputVector,
                                             vtkInformationVector *outputVector)
{
  //  Get the input and output data objects.
  //  Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //  Get the input
  vtkPolyData *polydata = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  /********************************************
   * Create a SurfaceMesh from the input mesh *
   ********************************************/
  SM sm;
  VPMap vpmap = get(CGAL::vertex_point, sm);
  
  //  Get nb of points and cells
  vtkIdType nb_points = polydata->GetNumberOfPoints();
  vtkIdType nb_cells = polydata->GetNumberOfCells();
  
  // Extract points
  std::vector<vertex_descriptor> vertex_map(nb_points);
  for (vtkIdType i=0; i<nb_points; ++i)
  {
    double coords[3];
    polydata->GetPoint(i, coords);
    vertex_descriptor v = add_vertex(sm);
    put(vpmap, v, K::Point_3(coords[0], coords[1], coords[2]));
    vertex_map[i] = v;
  }
  
  // Extract cells
  for (vtkIdType i = 0; i<nb_cells; ++i)
  {
    vtkCell* cell_ptr = polydata->GetCell(i);
    vtkIdType nb_vertices = cell_ptr->GetNumberOfPoints();
    std::vector<vertex_descriptor> vr(nb_vertices);
    for (vtkIdType k=0; k<nb_vertices; ++k)
      vr[k] = vertex_map[cell_ptr->GetPointId(k)];
    CGAL::Euler::add_face(vr, sm);
  }
  
  std::vector<vertex_descriptor> isolated_vertices;
  for(SM::vertex_iterator vit = sm.vertices_begin();
      vit != sm.vertices_end();
      ++vit)
  {
    if(sm.is_isolated(*vit))
      isolated_vertices.push_back(*vit);
  }
  
  for (std::size_t i=0; i < isolated_vertices.size(); ++i)
    sm.remove_vertex(isolated_vertices[i]);

  if(!is_triangle_mesh(sm))
  {
    vtkErrorMacro("The input mesh must be triangulated ");
    return 0;
  }
  
  /*****************************
   * Apply Isotropic remeshing *
   *****************************/
  PMP::isotropic_remeshing(sm.faces(),
                           Length,
                           sm,
                           PMP::parameters::number_of_iterations(MainIterations));
  
  /**********************************
   * Pass the SM data to the output *
   **********************************/
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkNew<vtkPoints> const vtk_points;
  vtkNew<vtkCellArray> const vtk_cells;
  vtk_points->Allocate(sm.number_of_vertices());
  vtk_cells->Allocate(sm.number_of_faces());
  std::vector<vtkIdType> Vids(num_vertices(sm));
  vtkIdType inum = 0;
  
  for(vertex_descriptor v : vertices(sm))
  {
    const K::Point_3& p = get(vpmap, v);
    vtk_points->InsertNextPoint(CGAL::to_double(p.x()),
                                CGAL::to_double(p.y()),
                                CGAL::to_double(p.z()));
    Vids[v] = inum++;
  }
  
  for(face_descriptor f : faces(sm))
  {
    vtkNew<vtkIdList> cell;
    for(halfedge_descriptor h : halfedges_around_face(halfedge(f, sm), sm))
      cell->InsertNextId(Vids[target(h, sm)]);

    vtk_cells->InsertNextCell(cell);
  }
  
  output->SetPoints(vtk_points);
  output->SetPolys(vtk_cells);
  output->Squeeze();

  return 1;
}

// ----------------------------------------------------------------------------
void vtkIsotropicRemeshingFilter::PrintSelf(std::ostream& os, 
                                            vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os<<"Length        : "<<Length        <<std::endl;
  os<<"LengthInfo    : "<<LengthInfo    <<std::endl;
  os<<"MainIterations: "<<MainIterations<<std::endl;
}

// ------------------------------------------------------------------------------
int vtkIsotropicRemeshingFilter::FillInputPortInformation(int vtkNotUsed(port),
                                                          vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

// ------------------------------------------------------------------------------
int vtkIsotropicRemeshingFilter::FillOutputPortInformation(int,
                                                           vtkInformation *info)
{
  // Always returns a vtkPolyData
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

// ------------------------------------------------------------------------------
int vtkIsotropicRemeshingFilter::RequestInformation(vtkInformation *,
                                                    vtkInformationVector ** inputVector,
                                                    vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  
  // Sets the bounds of the output.
  outInfo->Set(vtkDataObject::BOUNDING_BOX(),
               inInfo->Get(vtkDataObject::BOUNDING_BOX()),
               6);

  vtkPolyData *input= vtkPolyData::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  // Computes the initial target length:
  double * bounds = input->GetBounds();
  double diagonal = std::sqrt((bounds[0]-bounds[1]) * (bounds[0]-bounds[1]) +
                              (bounds[2]-bounds[3]) * (bounds[2]-bounds[3]) +
                              (bounds[4]-bounds[5]) * (bounds[4]-bounds[5]));
  SetLengthInfo(0.01*diagonal);
  
    return 1;
}
