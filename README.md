### Making a CGAL Plugin for Paraview


To make a plugin for Paraview, you need a developer version of Paraview, 
as a plugin only works with the same version it was built with.
To show you how, let us take an example: an Isotropic Remeshing plugin for Paraview. 

**ATTENTION**: This document is only valid until Paraview 5.6. The next version is currently (02/2019) under development, compatibility has been broken and documentation is not yet available. 


### Building Paraview from Sources

Go to [Paraview's GitLab page](https://gitlab.kitware.com/paraview/paraview) and get Paraview either by downloading the version you want, or by cloning the repository.
It is recommended to take a release version for safety and stability.
Once it is done, if you chose to clone the repository, you need to update the submodules. To do so,  go to the source dir, and do 

```
git submodule update --init --recursive
 ```
Once you have everything, create a build directory, and type this inside of it:
```
cmake <path-to-src>
make
```
It takes a little while, but in the end it should produce a usable version of Paraview.

### Writing a Plugin

Now that you have Paraview, assuming you already have a version of CGAL, you are ready to start coding.

The first thing to understand here is that VTK is based on a pipeline. 
As a quick overview, let's say that a SOURCE creates a VTK object, which can serve as input of a FILTER, that performs its algorithm when its `Update()` function is called. 


A Paraview plugin is the association between a filter of this pipeline and an XML file that is used to create a UI element.
Most of the time, we want to interact with a VTK data structure. In this example, our input is a `vtkPolyData` or a `vtkUnstructeredGrid`, so we choose to derive from a `vtkGeometryFilter`, 
that allows to convert a `vtkUnstructuredgrid` into a `vtkPolyData`.

A filter must implement the following methods to interact with the pipeline: 

- `RequestDataObject()`
   This function is used to create the output object. In this example, this is where the `vtkPolyData` that will hold the remeshed data is constructed.
   In this example, it is not needed, because the `vtkGeometryFilter` already implements it. It always returns an empty `vtkPolyData`.
- `RequestInformation()`
   In this function, we compute everything we can without heavy computation, like the bounding box for example. This is called before `Update()` and provides info that might be needed from the input. In this example, this is where we compute the default target edge length.
    Note: In some algorithms, this function is replaced by `ExecuteInformation()`. Before VTK 5, `Executeinformation()` was the standard, and some algorithm used for the transition kept it this way, so check it in the documentation of the base class you choose before you write that function, as it might never be called.
- `RequestData()`
   This is where the algorithm is performed. It gets the input, applies the algorithm and fills the output.
- `FillInputPortInformation()`
   This is where we specify what type of object we want as input.
- `FillOutputPortInformation()`
   This is where we specify what type of object will be outputted.



#### CMakeLists

On the cmake side, the most important thing is the `ADD_PARAVIEW_PLUGIN` macro. It deals with most of the VTK parts, provided the sources and the xml file for your plugin.

#### Header

#### Sources


#### XML
https://github.com/CGAL/cgal-paraview-plugins/blob/b80ba757641eef163ba4127754eb5799bfa2df72/IsotropicRemeshingFilter.xml#L1-L60

### Building the Plugin

To build a plugin, use CMake and provide the location of the CGAL library as CGAL_DIR.
Example : 
```
 cd /path/to/build/dir
 cmake -DCGAL_DIR=~/CGAL/releases/CGAL-4.9/cmake/platforms/release  -DParaView_DIR=<PATH_TO_Paraview_DIR>  path/to/dir/containing/plugin's/CMakeLists.txt
 make
```

Until Paraview 5.6.0 included, the Paraview_DIR is simply the build dir, where ParaViewConfig.cmake can be found.

Be careful to use the same version of Qt for compiling the plugin than the one used to compile paraview.

### Loading the Plugin in Paraview

Launch Paraview, go to Tools->Manage Plugins and click on Load New. Select the lib file of your plugin in the list and click Close. The plugin should appear in the Filter List.
