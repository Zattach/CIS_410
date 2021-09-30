/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
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
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)*idx[0];
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
//      vertexId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int vertexId, const int *dims)
{
    // 3D
    // idx[0] = vertexId%dim[0];
    // idx[1] = (vertexId/dims[0])%dims[1];
    // idx[2] = vertexId/(dims[0]*dims[1]);

    // 2D
    idx[0] = vertexId%dims[0];
    idx[1] = vertexId/dims[0];
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

// ****************************************************************************
//  Function: BoundingBoxForCell
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     cellId: a cellIndex (I.e., between 0 and GetNumberOfCells(dims))
//     bbox (output): the bounding box of cellId.  Format should be
//                     bbox[0]: the minimum X value in cellId.
//                     bbox[1]: the maximum X value in cellId.
//                     bbox[2]: the minimum Y value in cellId.
//                     bbox[3]: the maximum Y value in cellId.
//
//  Returns:  None (argument bbox is output)
//
// ****************************************************************************

void
BoundingBoxForCell(const float *X, const float *Y, const int *dims,
                   int cellId, float *bbox)
{
    int numberOfCells = GetNumberOfCells(dims);

    if(cellId >= 0 && cellId < numberOfCells){
        int logicalCellId[2];
        GetLogicalCellIndex(logicalCellId, cellId, dims);

        bbox[0] = X[logicalCellId[0]];
        bbox[1] = X[logicalCellId[0] + 1];
        bbox[2] = Y[logicalCellId[1]];
        bbox[3] = Y[logicalCellId[1] + 1];
    } else {
        bbox[0] = -100;
        bbox[1] = 100;
        bbox[2] = -100;
        bbox[3] = 100;
    }
}

// ****************************************************************************
//  Function: CountNumberOfStraddlingCells
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//  Returns:  the number of cells that straddle 0, i.e., the number of cells
//            that contains points who have F>0 and also have points with F<0.
//
// ****************************************************************************

int
CountNumberOfStraddlingCells(const float *X, const float *Y, const int *dims,
                             const float *F)
{
    int numberOfStraddlingCells = 0;
    int numberOfCells = GetNumberOfCells(dims);
    int logicalCellId[2];
    int vertexId[4];

    for(int i=0; i < numberOfCells; i++){
        GetLogicalCellIndex(logicalCellId, i, dims);

        vertexId[0] = GetPointIndex(logicalCellId, dims);

        logicalCellId[0]++;
        vertexId[1] = GetPointIndex(logicalCellId, dims);
        
        logicalCellId[0]--; logicalCellId[1]++;
        vertexId[2] = GetPointIndex(logicalCellId, dims);
        
        logicalCellId[0]++;
        vertexId[3] = GetPointIndex(logicalCellId, dims);

        if(!(F[vertexId[0]] > 0 && F[vertexId[1]] > 0 && F[vertexId[2]] > 0 && F[vertexId[3]] > 0) 
            && !(F[vertexId[0]] < 0 && F[vertexId[1]] < 0 && F[vertexId[2]] < 0 && F[vertexId[3]] < 0)){
            numberOfStraddlingCells++;
        }
    }
    return numberOfStraddlingCells;
}

int main()
{
    int  i;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj2_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int numCells = CountNumberOfStraddlingCells(X, Y, dims, F);
    cerr << "The number of cells straddling zero is " << numCells << endl;

    float bbox[4];
    const int ncells = 5;
    int cellIds[ncells] = { 0, 50, 678, 1000, 1200 };
    for (i = 0 ; i < ncells ; i++)
    {
        BoundingBoxForCell(X, Y, dims, cellIds[i], bbox);
        cerr << "The bounding box for cell " << cellIds[i] << " is " 
             << bbox[0] << "->" << bbox[1] << ", " << bbox[2] << "->" << bbox[3]
             << endl;
    }

    const int npts = 10;
    float pt[npts][3] = 
         {
            {1.01119, 0.122062, 0},
            {0.862376, 1.33839, 0},
            {0.155026, 0.126123, 0},
            {0.69736, 0.0653565, 0},
            {0.2, 0.274117, 0},
            {0.893699, 1.04111, 0},
            {0.608791, -0.0533753, 0},
            {1.00543, 0.138024, 0},
            {0.384128, -0.0768977, 0},
            {0.666757, 0.60259, 0},
         };

    

    for (i = 0 ; i < npts ; i++)
    {
        float f = EvaluateFieldAtLocation(pt[i], dims, X, Y, F);
        cerr << "Evaluated field at (" << pt[i][0] <<"," << pt[i][1] << ") as "
             << f << endl;
    }
}