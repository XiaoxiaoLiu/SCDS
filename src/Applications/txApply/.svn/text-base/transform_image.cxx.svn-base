#include <assert.h>

#include <iostream>
#include <limits>
#include <string>

#include "Array3D.h"
#include "Array3DIO.h"
#include "Array3DUtils.h"
#include "Debugging.h"
#include "HField3DIO.h"
#include "HField3DUtils.h"
#include "Image.h"
#include "ImageIO.h"
#include "ImageUtils.h"

typedef float VoxelType;
float MAX_VOX = std::numeric_limits<float>::max();

const char* defaultProgramName = "transform_image";
const char* programName = defaultProgramName;

void loadImage(const char *fileName, Image<VoxelType>& image)
{
  ImageIO imageLoader;
  try
    {
      imageLoader.LoadThisImage(fileName, image);
    }
  catch (...)
    {
      std::cerr << "Error loading file: " << fileName << std::endl;
      exit(0);
    }
}

void printUsage(const std::string& message)
{
    std::cerr 
        << message
        << "\nUsage: " << programName << " followed by:" 
        << "\n  -i InputImageName"
        << "\n  -o Output3DImagename"
        << "\n  -t TemplateName (image to use as template for origin and spacing)"
        << "\n  -h HFieldName (deformation field)"
        << "\n  -m MatrixName (plunc format)"
        << "\n  -x SliceNum SliceAxis PgmImageName (see below)"
        << "\n  -w Min Max (intensity window)"
        << "\n"
        << "\nThe -x option outputs a pgm image to file named PgmImageName of the"
        << "\nslice with index SliceNum, perpendicular to the axis determined by"
        << "\nSliceAxis.  For SliceAxis, 0 means sagittal, 1 means coronal, and 2"
        << "\nmeans axial.  The -x option can appear multiple times to output"
        << "\ndifferent slices."
        << "\n"
        << "\nThe only required option is -i.  There can be multiple -h"
        << "\nand -m options; if none, the identity transform is assumed."
        << "\nIf there is no -t option, the original image origin and "
        << "\nspacing is used."
        << std::endl;
}

void
WritePgmSlice(const std::string& filename, const Image<VoxelType>& image,
              int slice, int axis, float intensMin = 0, float intensMax = 0)
{

    Vector3D<unsigned int> end = image.getSize();
    end[axis] = slice + 1;
    Vector3D<unsigned int> begin(0, 0, 0);
    begin[axis] = slice;

    std::cout << begin << ";  " << end << std::endl;

    int numPixels = (end[0] - begin[0]) * (end[1] - begin[1]) * (end[2] - begin[2]);
    std::cout << OH(numPixels) << std::endl;
    unsigned char* pixels = new unsigned char[numPixels];

    unsigned int k, j, i;

    // Do this only if you need to find the max and min over that slice
    if (intensMin == intensMax) {
        intensMin = MAX_VOX;
        intensMax = -MAX_VOX;
        for (k = begin[2]; k < end[2]; ++k) {
            for (j = begin[1]; j < end[1]; ++j) {
                for (i = begin[0]; i < end[0]; ++i) {
                    float pixVal = image(i, j, k);
                    if (pixVal < intensMin) intensMin = pixVal;
                    if (pixVal > intensMax) intensMax = pixVal;
                }
            }
        }
    }

    unsigned char* pix = pixels;
    for (k = begin[2]; k < end[2]; ++k) {
        for (j = begin[1]; j < end[1]; ++j) {
            for (i = begin[0]; i < end[0]; ++i) {
                float pixVal = (image(i, j, k) - intensMin) / 
                               (intensMax - intensMin) * 255.99f;
                if (pixVal < 0) pixVal = 0;
                if (pixVal > 255) pixVal = 255;
                *pix++ = (unsigned char) pixVal;
            }
        }
    }

    // The begin[i] terms should be redundant here, but we might put in
    // ROI stuff.
    int width, height;
    if (axis == 0) { width = end[1] - begin[1]; height = end[2] - begin[2]; } 
    if (axis == 1) { width = end[0] - begin[0]; height = end[2] - begin[2]; } 
    if (axis == 2) { width = end[0] - begin[0]; height = end[1] - begin[1]; } 

    FILE* fp = fopen(filename.c_str(), "wb");
    fprintf(fp, "P5\n%d %d\n255\n", width, height);
    fwrite(pixels, 1, numPixels, fp);
    fclose(fp);
}


int main(int argc, char **argv)
{

    std::string outputFilename = "";
    std::string imageFilename  = "";
    std::string templateFilename  = "";

    Vector3D<unsigned int> size;
    Vector3D<double> origin;
    Vector3D<double> spacing;

    std::vector<char*> transformTypes;
    std::vector<char*> transformFilenames;

    std::vector< std::string > sliceFilenames;
    std::vector<int> sliceAxes;
    std::vector<int> slicesToSave;
    float intensMin = 0;
    float intensMax = 0;

    programName = argv[0];

    argv++; argc--;
    while (argc > 0) {
        std::string thisArg(argv[0]);
        if (thisArg == "-o") {
            // outfile
            if (argc < 2) { printUsage("infile"); return -1; }
            outputFilename = argv[1];
            argv += 2; argc -= 2;
	} else if (std::string(argv[0]) == std::string("-i")) {
            // image
            if (argc < 2) { printUsage("outfile"); return -1; }
            imageFilename = argv[1];
            argv += 2; argc -= 2;
	} else if (std::string(argv[0]) == std::string("-t")) {
            // template for origin & spacing
            if (argc < 2) { printUsage("template"); return -1; }
            templateFilename = argv[1];
            argv += 2; argc -= 2;
	} else if (std::string(argv[0]) == std::string("-h") 
                   || std::string(argv[0]) == std::string("-m") ) {
            // hfield or matrix
            if (argc < 2) { printUsage("transform"); return -1; }
            transformTypes.push_back(argv[0]);
            transformFilenames.push_back(argv[1]);
            argv += 2; argc -= 2;
	} else if (std::string(argv[0]) == std::string("-x")) {
            // save slice
            if (argc < 4) { printUsage("slice"); return  -1; }
            slicesToSave.push_back( atoi(argv[1]) );
            sliceAxes.push_back( atoi(argv[2]) );
            sliceFilenames.push_back( std::string(argv[3]) );
            argv += 4; argc -= 4;
	} else if (std::string(argv[0]) == std::string("-w")) {
            // set intensity window
            if (argc < 3) { printUsage("intensity window"); return  -1; }
            intensMin = atof(argv[1]);
            intensMax = atof(argv[2]);
            argv += 3; argc -= 3;
        } else { printUsage(argv[0]); return 0; }
    }

    if (imageFilename == "") { printUsage("No input image."); return 0; }
    assert(transformTypes.size() == transformFilenames.size());

    Image<float> image;
    std::cout << "Loading image: " << imageFilename << "...";
    loadImage(imageFilename.c_str(), image);
    std::cout << "DONE" << std::endl;
    std::cout << "  Dimensions: " << image.getSize() << std::endl;
    std::cout << "  Origin    : " << image.getOrigin() << std::endl;
    std::cout << "  Spacing   : " << image.getSpacing() << std::endl;
    Image<float> def;

    if (templateFilename != "") {
        ImageIO headerLoader;
        headerLoader.LoadHeader(templateFilename, size, origin, spacing);
    } else {
        size = image.getSize();
        origin = image.getOrigin();
        spacing = image.getSpacing();
    }

    std::cout << "\nspacing = " << spacing << std::endl;

    // Apply transformations
    for (unsigned int i = 0; i < transformTypes.size(); ++i) {

        if (std::string(transformTypes[i]) == "-h") { // Hfield

            Array3D<Vector3D<float> > hField;
            Vector3D<unsigned int> hSize;
            Vector3D<double> hOrigin, hSpacing;
            const char* hff = transformFilenames[i];
            std::cout << "Loading hField: " << hff << "..." << std::flush;

            HField3DIO::readMETAHeader(hSize, hOrigin, hSpacing, hff);
            HField3DIO::readMETA(hField, hff);

            std::cout << "Done." 
                      << "\n  Dimensions: " << hSize 
                      << "\n  Origin    : " << hOrigin 
                      << "\n  Spacing   : " << hSpacing << std::endl;

            HField3DUtils::apply(image, hField, def, hOrigin, hSpacing);
            ImageUtils::resampleWithTransparency(def, image);

        } else if (std::string(transformTypes[i]) == "-m") { // matrix
            AffineTransform3D<float> affine;
            affine.readPLUNCStyle(transformFilenames[i]);
            ImageUtils::applyAffine(image, origin, spacing, size, affine);
        }
    }

    assert(slicesToSave.size() == sliceAxes.size() &&
           slicesToSave.size() == sliceFilenames.size());
    for (unsigned int i = 0; i < slicesToSave.size(); ++i) {
        WritePgmSlice(sliceFilenames[i], image, slicesToSave[i], sliceAxes[i],
                      intensMin, intensMax);
    }

    if (outputFilename != "" ) {
        ImageIO imio;
        imio.SaveImage(outputFilename, image);
    }
    return 0;
}
  
  
  
