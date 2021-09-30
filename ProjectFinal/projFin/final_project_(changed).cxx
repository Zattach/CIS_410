#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <math.h>

#define PI 3.14159265


struct Camera
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
};


struct Ray
{
    double origin[3];
    double direction[3];
};


struct TransferFunction
{
    double          min;
    double          max;
    int             numBins;
    unsigned char  *colors;  // size is 3*numBins
    double         *opacities; // size is numBins


    // Finds which bin
    int GetBin(double fieldValue){
        if(fieldValue < min || fieldValue > max){
            // cout << "Out of range: no valid bin" << endl;    // output message
            return -1;
        }

        // Locate bin using log(n) time
        int low = 0;
        int high = numBins - 1;
        int bin;
        bool found = false;
        while(!found){
            bin = ((high + low) / 2);
            double lowFieldVal = min + (((max - min) / numBins) * bin);
            double highFieldVal = min + (((max - min) / numBins) * (bin + 1));
            if (fieldValue >= lowFieldVal && fieldValue < highFieldVal){
                found = true;
            } else if(fieldValue < lowFieldVal){
                high = bin - 1;
            } else {
                low = bin + 1;
            }
        }

        // cout << "Mapped to bin " << bin << endl;     // output message
        return bin;
    }


    // Applies the transfer function to the given fieldValue
    void ApplyTransferFunction(double value, unsigned char *RGB, double &opacity)
    {
        int bin = GetBin(value);

        // Manually sets all values to 0 iff the fieldValue is out of the range
        if(bin == -1){
            RGB[0] = 0;
            RGB[1] = 0;
            RGB[2] = 0;
            opacity = 0;
        } else {
            RGB[0] = colors[3*bin+0];
            RGB[1] = colors[3*bin+1];
            RGB[2] = colors[3*bin+2];
            opacity = opacities[bin];
        }
    }
};


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
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    // return idx[1]*dims[0]+idx[0];
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

double
EvaluateFieldAtLocation(const double *pt, const int *dims, 
                        const float *X, const float *Y, const float *Z, const float *F)
{
    // Check if pt[] is a valid point
    float tempPt[3] = {(float) pt[0], (float) pt[1], (float) pt[2]};
    if(tempPt[0] < X[0] || tempPt[1] < Y[0] || tempPt[2] < Z[0] || 
        tempPt[0] > X[dims[0] - 1] || tempPt[1] > Y[dims[1] - 1] || tempPt[2] > Z[dims[2] - 1]){
        return 0;
    }

    int firstVertex[3];      // Bottom left vertex of cell containing interpolated point

    // Locate X value for the bottom left vertex of the cell containing pt[]
    int low = 0;
    int high = dims[0] - 1;
    bool foundX = false;
    while(!foundX){
        int mid = (low + high) / 2;
        if (tempPt[0] >= X[mid] && tempPt[0] <= X[mid + 1]){
            foundX = true;
            firstVertex[0] = mid;
        } else if(tempPt[0] < X[mid]){
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
        if (tempPt[1] >= Y[mid] && tempPt[1] < Y[mid + 1]){
            foundY = true;
            firstVertex[1] = mid ;
        } else if(tempPt[1] < Y[mid]){
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }

    // Locate Z value for the bottom left vertex of the cell containing pt[]
    low = 0;
    high = dims[1] - 1;
    bool foundZ = false;
    while(!foundZ){
        int mid = (low + high) / 2;
        if (tempPt[2] >= Z[mid] && tempPt[2] < Z[mid + 1]){
            foundZ = true;
            firstVertex[2] = mid ;
        } else if(tempPt[2] < Z[mid]){
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }


    // Find point index values for each vertex of the cell containing pt[]
    int vertexId[8];
    vertexId[0] = GetPointIndex(firstVertex, dims);
    firstVertex[0]++;
    vertexId[1] = GetPointIndex(firstVertex, dims);
    firstVertex[0]--;    firstVertex[1]++;
    vertexId[2] = GetPointIndex(firstVertex, dims);
    firstVertex[0]++;
    vertexId[3] = GetPointIndex(firstVertex, dims);
    firstVertex[0]--;    firstVertex[1]--;      firstVertex[2]++;
    vertexId[4] = GetPointIndex(firstVertex, dims);
    firstVertex[0]++;
    vertexId[5] = GetPointIndex(firstVertex, dims);
    firstVertex[0]--;    firstVertex[1]++;
    vertexId[6] = GetPointIndex(firstVertex, dims);
    firstVertex[0]++;
    vertexId[7] = GetPointIndex(firstVertex, dims);
    firstVertex[0]--;    firstVertex[1]--;      firstVertex[2]--;

    // Calculate necessary t-values for both X and Y
    double tValueX = (double) (tempPt[0] - X[firstVertex[0]]) / (X[firstVertex[0] + 1] - X[firstVertex[0]]);
    double tValueY = (double) (tempPt[1] - Y[firstVertex[1]]) / (Y[firstVertex[1] + 1] - Y[firstVertex[1]]);
    double tValueZ = (double) (tempPt[2] - Z[firstVertex[2]]) / (Z[firstVertex[2] + 1] - Z[firstVertex[2]]);

    // Calculate interpolated F-values for temporary points
    double interpolatedLowerValue = (double) (F[vertexId[0]] + (tValueX * (F[vertexId[1]] - F[vertexId[0]])));
    double interpolatedHigherValue = (double) (F[vertexId[2]] + (tValueX * (F[vertexId[3]] - F[vertexId[2]])));

    // Calculate interpolated F-values for front vertex
    double interpolatedFrontPtValue = (double) interpolatedLowerValue + (tValueY * (interpolatedHigherValue - interpolatedLowerValue));

    // Calculate interpolated F-values for Z += 1
    interpolatedLowerValue = (double) (F[vertexId[4]] + (tValueX * (F[vertexId[5]] - F[vertexId[4]])));
    interpolatedHigherValue = (double) (F[vertexId[6]] + (tValueX * (F[vertexId[7]] - F[vertexId[6]])));

    // Calculate interpolated F-values for back vertex
    double interpolatedBackPtValue = interpolatedLowerValue + (((double) tValueY) * (interpolatedHigherValue - interpolatedLowerValue));

    // Calculates final interpolated F-value
    double interpolatedPtValue = interpolatedFrontPtValue + (((double) tValueZ) * (interpolatedBackPtValue - interpolatedFrontPtValue));
    return interpolatedPtValue;
}


void WriteImage(vtkImageData *img, const char *filename){
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData * NewImage(int width, int height){
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}


TransferFunction SetupTransferFunction(void){
    int  i;

    TransferFunction rv;
    rv.min = 10;
    rv.max = 15;
    rv.numBins = 256;
    rv.colors = new unsigned char[3*256];
    rv.opacities = new double[256];
    unsigned char charOpacity[256] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };

    for (i = 0 ; i < 256 ; i++)
        rv.opacities[i] = charOpacity[i]/255.0;
    const int numControlPoints = 8;
    unsigned char controlPointColors[numControlPoints*3] = { 
           71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0, 
           255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76 
       };
    double controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
    for (i = 0 ; i < numControlPoints-1 ; i++)
    {
        int start = controlPointPositions[i]*rv.numBins;
        int end   = controlPointPositions[i+1]*rv.numBins+1;
cerr << "Working on " << i << "/" << i+1 << ", with range " << start << "/" << end << endl;
        if (end >= rv.numBins)
            end = rv.numBins-1;
        for (int j = start ; j <= end ; j++)
        {
            double proportion = (j/(rv.numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
            if (proportion < 0 || proportion > 1.)
                continue;
            for (int k = 0 ; k < 3 ; k++)
                rv.colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k])
                                 + controlPointColors[3*i+k];
        }
    }    

    return rv;
}


Camera
SetupCamera(void)
{
    Camera rv;
    rv.focus[0] = 0;
    rv.focus[1] = 0;
    rv.focus[2] = 0;
    rv.up[0] = 0;
    rv.up[1] = -1;
    rv.up[2] = 0;
    rv.angle = 30;
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position[0] = -8.25e+7;
    rv.position[1] = -3.45e+7;
    rv.position[2] = 3.35e+7;

    return rv;
}

void crossProduct(double vectA[], double vectB[], double crossP[]){
    crossP[0] = vectA[1] * vectB[2] - vectA[2] * vectB[1];
    crossP[1] = vectA[2] * vectB[0] - vectA[0] * vectB[2];
    crossP[2] = vectA[0] * vectB[1] - vectA[1] * vectB[0];
}

Ray getRay(Camera cam, int nx, int ny, int px, int py){
    Ray r;
    r.origin[0] = cam.position[0];
    r.origin[1] = cam.position[1];
    r.origin[2] = cam.position[2];

    double look[3];
    for(int i = 0; i < 3; i++)
        look[i] = cam.focus[i] - cam.position[i];

    double temp[3];
    double length;
    crossProduct(look, cam.up, temp);
    length = sqrt((temp[0] * temp[0]) + (temp[1] * temp[1]) + (temp[2] * temp[2]));
    double R_u[3] = {temp[0] / length, temp[1] / length, temp[2] / length};

    crossProduct(look, R_u, temp);
    length = sqrt((temp[0] * temp[0]) + (temp[1] * temp[1]) + (temp[2] * temp[2]));
    double R_v[3] = {temp[0] / length, temp[1] / length, temp[2] / length};

    double tempScalar = (2 * tan((cam.angle * PI / 180) / 2)) / nx;
    double R_x[3] = {tempScalar * R_u[0], tempScalar * R_u[1], tempScalar * R_u[2]};

    tempScalar = (2 * tan((cam.angle * PI / 180) / 2)) / ny;
    double R_y[3] = {tempScalar * R_v[0], tempScalar * R_v[1], tempScalar * R_v[2]};

    length = sqrt((look[0] * look[0]) + (look[1] * look[1]) + (look[2] * look[2]));
    for(int i = 0; i < 3; i++)
        r.direction[i] = (look[i] / length) + (((2.0 * px + 1 - nx) / 2.0) * R_x[i]) + (((2.0 * py + 1 - ny) / 2.0) * R_y[i]);

    return r;
}


int main()
{
    int nx = 1000;
    int ny = 1000;
    int sampleSize = 1024;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("astro512.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    vtkImageData *image = NewImage(nx, ny);
    unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);

    for (int i = 0 ; i < 3*nx*ny ; i++)
        buffer[i] = 0;

    // Declare necessary variables
    Camera cam = SetupCamera();
    TransferFunction tf = SetupTransferFunction();
    Ray r_d;
    float stepSize = (cam.far - cam.near) / (sampleSize - 1);
    double pt[3];

    // Loop through every pixel to determine the RGB value based on the ray being cast through the volume
    for(int i = 0 ; i < nx ; i++){
        for(int j = 0 ; j < ny ; j++){
            cout << "Solving point " << i << ", " << j << endl;     // Output message to see progress (not frozen)

            // Get the origin and direction of the ray for pixel i, j
            r_d = getRay(cam, nx, ny, i, j);

            // Set all necessary values to 0
            double RGBRun[3] = {0, 0, 0};
            unsigned char RGBCurr[3] = {0, 0, 0};
            double opacityPrev = 0.0;
            double opacityCurr = 0.0;
            double opacityRun = 0.0;
            for(int sample = 0; sample < sampleSize; sample++){

                // Find the location for the current sample point
                for(int sample_dim = 0; sample_dim < 3; sample_dim++)
                    pt[sample_dim] = r_d.origin[sample_dim] + (r_d.direction[sample_dim] * cam.near) + (r_d.direction[sample_dim] * stepSize * sample);

                // Evaluate the field value of the sample point and apply it to the transfer function
                double fieldValue = EvaluateFieldAtLocation(pt, dims, X, Y, Z, F);
                tf.ApplyTransferFunction(fieldValue, RGBCurr, opacityCurr);

                // Correct the opacity
                opacityCurr = 1.0 - pow(1.0 - opacityCurr, 500.0 / sampleSize);

                // Color composite over the next sample of the ray
                opacityPrev = opacityRun;
                opacityRun = opacityRun + ((1 - opacityRun) * opacityCurr);
                if(opacityRun == 0){
                    RGBRun[0] = 0;
                    RGBRun[1] = 0;
                    RGBRun[2] = 0;
                } else {
                    RGBRun[0] = ((opacityPrev * RGBRun[0]) + ((1 - opacityPrev) * opacityCurr * (RGBCurr[0] / 255.0))) / opacityRun;
                    RGBRun[1] = ((opacityPrev * RGBRun[1]) + ((1 - opacityPrev) * opacityCurr * (RGBCurr[1] / 255.0))) / opacityRun;
                    RGBRun[2] = ((opacityPrev * RGBRun[2]) + ((1 - opacityPrev) * opacityCurr * (RGBCurr[2] / 255.0))) / opacityRun;
                }
            }

            // Color composite over a black background
            RGBCurr[0] = 0;
            RGBCurr[1] = 0;
            RGBCurr[2] = 0;
            opacityCurr = 1;

            opacityPrev = opacityRun;
            opacityRun = opacityRun + ((1 - opacityRun) * opacityCurr);
            if(opacityRun == 0){
                RGBRun[0] = 0;
                RGBRun[1] = 0;
                RGBRun[2] = 0;
            } else {
                RGBRun[0] = ((opacityPrev * RGBRun[0]) + ((1 - opacityPrev) * opacityCurr * (RGBCurr[0] / 255.0))) / opacityRun;
                RGBRun[1] = ((opacityPrev * RGBRun[1]) + ((1 - opacityPrev) * opacityCurr * (RGBCurr[1] / 255.0))) / opacityRun;
                RGBRun[2] = ((opacityPrev * RGBRun[2]) + ((1 - opacityPrev) * opacityCurr * (RGBCurr[2] / 255.0))) / opacityRun;
            }

            // Save RGB values to the buffer
            int offset = 3*(j*nx+i);
            buffer[offset] = (unsigned char) (RGBRun[0] * 255.0);
            buffer[offset + 1] = (unsigned char) (RGBRun[1] * 255.0);
            buffer[offset + 2] = (unsigned char) (RGBRun[2] * 255.0);
        }
    }

    WriteImage(image, "my_astro");
}
