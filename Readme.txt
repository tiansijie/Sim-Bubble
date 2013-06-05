This readme file is intended to explain how to build the project. To use the authoring tool, please refer to Tutorial section in the Final Report.

The solution can currently only be built as release version as no debug version library file for CGAL is built in debug version.
There are two projects within solution. Only 'BubbleSimMaya' is required to build to generate mll file, 'Bubble' is the OpenGL test program.

In order for the 'BubbleSimMaya' project to build, in project property page, do the following:

In c/c++ -> General -> Additional Include Directories, 1. change D:\Program Files\Autodesk\Maya2013\include to local Maya include directory 2. change E:\Program Files\boost\boost_1_51 to boost root directory (boost 1.51 is used here) 3. change G:\Visual Studio 2010\Projects\Advanced CG\Bubble\gmp\include to the gmp include folder in solution directory 4. change G:\Visual Studio 2010\Projects\Advanced CG\Bubble\CGAL Install\include to CGAL include folder in solution directory


In Linker -> General ->Additional Library Directories, 1. change D:\Program Files\Autodesk\Maya2013\lib to local Maya lib directory
2. change E:\Program Files\boost\boost_1_51\lib to boost lib directory (boost 1.51 is used here) 3. change G:\Visual Studio 2010\Projects\Advanced CG\Bubble\gmp\lib to the gmp lib folder in solution directory 4. change G:\Visual Studio 2010\Projects\Advanced CG\Bubble\CGAL Install\lib to CGAL lib folder in solution directory


Files used in Tutorial is included in 'Auxiliary Files' directory