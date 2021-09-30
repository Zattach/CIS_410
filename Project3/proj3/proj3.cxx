#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float
EvaluateFieldAtLocation(const float *pt, const int *dims, 
                        const float *X, const float *Y, const float *F)
{
    // Check if pt[] is a valid point
    if(pt[0] < X[0] || pt[1] < Y[0] || pt[0] > X[dims[0] - 1] || pt[1] > Y[dims[1] - 1]){
        return 0;
    }

    int firstVertex[2];      // Bottom left vertex of cell containing interpolated point

    // Locate X value for the bottom left vertex of the cell containing pt[]
    int low = 0;
    int high = dims[0] - 1;
    bool foundX = false;
    while(!foundX){
        int mid = (low + high) / 2;
        if (pt[0] >= X[mid] && pt[0] < X[mid + 1]){
            foundX = true;
            firstVertex[0] = mid;
        } else if(pt[0] < X[mid]){
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }

    // Locate Y value for the bottom left vertex of the cell containing pt[]
    low = 0;
    high = dims[1] - 1;
    bool foundY = false;
    while(!foundY){
        int mid = (low + high) / 2;
        if (pt[1] >= Y[mid] && pt[1] < Y[mid + 1]){
            foundY = true;
            firstVertex[1] = mid ;
        } else if(pt[1] < Y[mid]){
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }

    // Find point index values for each vertex of the cell containgin pt[]
    int vertexId[4];
    vertexId[0] = GetPointIndex(firstVertex, dims);
    firstVertex[0]++;
    vertexId[1] = GetPointIndex(firstVertex, dims);
    firstVertex[0]--;    firstVertex[1]++;
    vertexId[2] = GetPointIndex(firstVertex, dims);
    firstVertex[0]++;
    vertexId[3] = GetPointIndex(firstVertex, dims);
    firstVertex[0]--;    firstVertex[1]--;

    // Calculate necessary t-values for both X and Y
    float tValueX = (pt[0] - X[firstVertex[0]]) / (X[firstVertex[0] + 1] - X[firstVertex[0]]);
    float tValueY = (pt[1] - Y[firstVertex[1]]) / (Y[firstVertex[1] + 1] - Y[firstVertex[1]]);

    // Calculate interpolated F-values for temporary points
    float interpolatedLowerValue = F[vertexId[0]] + (tValueX * (F[vertexId[1]] - F[vertexId[0]]));
    float interpolatedHigherValue = F[vertexId[2]] + (tValueX * (F[vertexId[3]] - F[vertexId[2]]));

    // Calculate interpolated F-values for final vertex
    float interpolatedPtValue = interpolatedLowerValue + (tValueY * (interpolatedHigherValue - interpolatedLowerValue));

    return interpolatedPtValue;
}


void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,128) 
//        F=1: (255,255,255) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************

void
ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
    RGB[0] = F * (255);
    RGB[1] = F * (255);
    RGB[2] = 128 + (F * (255 - 128));
}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using a divergent colormap
//
//     The divergent color map has:
//        F=0: (0,0,128) 
//        F=0.5: (255,255,255) 
//        F=1: (128, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
    if(F < 0.5){
        float tValue = 2 * F;
        RGB[0] = tValue * 255.0;
        RGB[1] = tValue * 255.0;
        RGB[2] = 128.0 + (tValue * (255.0 - 128.0));
    } else {
        float tValue = 2 * (F - .5);
        RGB[0] = 255.0 + (tValue * (128.0 - 255.0));
        RGB[1] = 255.0 + (tValue * -255.0);
        RGB[2] = 255.0 + (tValue * -255.0);
    }
}

// ****************************************************************************
//  Function: ApplyBHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using an HSV rainbow colormap
//
//     The rainbow colormap uses a saturation =1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees 
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyHSVColorMap(float F, unsigned char *RGB)
{
    float hue = F * 360.0;
    float saturation = 1.0;
    float value = 1.0;

    hue /= 60.f;

    // sector 0 to 5
    int initialValue = floor( hue );
    float remainder = hue - initialValue;
    
    // factorial part of h
    float p = value * ( 1.0 - saturation);
    float q = value * ( 1.0 - saturation * remainder );
    float t = value * ( 1.0 - saturation * ( 1.0 - remainder ) );

    switch(initialValue){
        case 0: RGB[0] = value * 255.0;
            RGB[1] = t * 255.0;
            RGB[2] = p * 255.0;
            break;
        case 1: RGB[0] = q * 255.0;
            RGB[1] = value * 255.0;
            RGB[2] = p * 255.0;
            break;
        case 2: RGB[0] = p * 255.0;
            RGB[1] = value * 255.0;
            RGB[2] = t * 255.0;
            break;
        case 3: RGB[0] = p * 255.0;
            RGB[1] = q * 255.0;
            RGB[2] = value * 255.0;
            break;
        case 4: RGB[0] = t * 255.0;
            RGB[1] = p * 255.0;
            RGB[2] = value * 255.0;
            break;
        default: RGB[0] = value * 255.0;
            RGB[1] = p * 255.0;
            RGB[2] = q * 255.0;
            break;
    }
}


int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int nx = 500;
    int ny = 500;

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3*nx*ny ; i++)
        for (j = 0 ; j < 3 ; j++)
            buffer[j][i] = 0;

    for (i = 0 ; i < nx ; i++)
        for (j = 0 ; j < ny ; j++)
        {
            // ITERATE OVER PIXELS
            float pt[2];
            pt[0] = ((18.0 * i) / (nx - 1)) - 9.0;
            pt[1] = ((18.0 * j) / (ny - 1)) - 9.0;
            float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);
            float normalizedF = (f - 1.2) / (5.02 - 1.2); // see step 5 re 1.2->5.02
            
            // I TAKE OVER HERE
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
            ApplyHSVColorMap(normalizedF, buffer[2]+offset);
        }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}
