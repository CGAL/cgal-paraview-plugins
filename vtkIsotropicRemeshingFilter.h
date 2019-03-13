#ifndef vtkIsotropicRemeshingFilter_h
#define vtkIsotropicRemeshingFilter_h
// Gives access to macros for communication with the UI
#include "vtkFiltersCoreModule.h" 
#include "vtkGeometryFilter.h"

// Inherit from the desired filter
class vtkIsotropicRemeshingFilter : public vtkGeometryFilter
{
public:
  // VTK requirements
  static vtkIsotropicRemeshingFilter* New();
  vtkTypeMacro(vtkIsotropicRemeshingFilter, vtkGeometryFilter);
  // Prints the values of the specific data
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Communicate with the UI
  vtkSetMacro(Length, double);
  vtkGetMacro(Length, double);
  vtkSetMacro(LengthInfo, double);
  vtkGetMacro(LengthInfo, double);
  vtkSetMacro(MainIterations, int);
  vtkGetMacro(MainIterations, int);

  // Pipeline functions:
  // Performs the isotropic remeshing algorithm and fills the output object here.
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *)override;
  // Specifies the type of the input objects
  int FillInputPortInformation(int, vtkInformation *info)override;
  // Specifies the type of the output object.
  int FillOutputPortInformation(int, vtkInformation *info)override;

protected:
  vtkIsotropicRemeshingFilter();
  ~vtkIsotropicRemeshingFilter(){}

  // Computes the bbox's diagonal length to set the default target edge length.
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  // Data set by the UI and used by the algorithm
  double Length;
  double LengthInfo;
  int MainIterations;

  // needed but not implemented
  vtkIsotropicRemeshingFilter(const vtkIsotropicRemeshingFilter&);
  void operator=(const vtkIsotropicRemeshingFilter&);
};
#endif
