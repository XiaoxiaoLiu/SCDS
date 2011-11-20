1.Get NuroLib
cvs -d :pserver:anonymous@demeter.ia.unc.edu:/cvsroot login (press Enter for password) 
cvs -d :pserver:anonymous@demeter.ia.unc.edu:/cvsroot co -P NeuroLib



2.install GenerateCLP package:
http://www.slicer.org/slicerWiki/index.php/Slicer3:Execution_Model_Documentation#Using_GenerateCLP_Outside_of_Slicer3

GenerateCLP can be built and used outside of the Slicer3 tree. These instructions show how to build GenerateCLP for your own projects without needing a full Slicer3 build or installation. (To see how to simply a build a Slicer3 module against a Slicer3 build or installation, see User-defined Modules Built Outside of Slicer3.
Check out the required directories from the Slicer3 repository.
svn co http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/tclap
svn co http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/ModuleDescriptionParser
svn co http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/GenerateCLP
Run cmake on tclap
Run cmake on ModuleDescriptionParser
make ModuleDescriptionParser
Run cmake on GenerateCLP
make GenerateCLP


3. My exprience in building postprocessingTPS in windows
(VS 2010 express)
 To build PostprocessingTPs, I have to build "static" library for all the depending libraries, including ITK, tclap, ModuleDesscriptionParse, GenerateCLP.
 Somehow if I use the default "shared" library, GenerateCLP can not generate ***CLP.h correctly. 
 
