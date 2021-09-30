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
#include <vtkPNGWriter.h>

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

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>


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


class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj5.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!

    // Populate numSegments array
    int numSegments[16];
    numSegments[0] = 0;
    numSegments[1] = numSegments[2] = numSegments[3] = numSegments[4] = numSegments[5] = 1;
    numSegments[6] = 2;
    numSegments[7] = numSegments[8] = 1;
    numSegments[9] = 2;
    numSegments[10] = numSegments[11] = numSegments[12] = numSegments[13] = numSegments[14] = 1;
    numSegments[15] = 0;

    // Populate lup table with data corresponding to cell templates
    int lup[16][4];
    lup[0][0] = lup[0][1] = lup[0][2] = lup[0][3] = -1;
    lup[1][0] = 0; lup[1][1] = 3; lup[1][2] = lup[1][3] = -1;
    lup[2][0] = 0; lup[2][1] = 1; lup[2][2] = lup[2][3] = -1;
    lup[3][0] = 1; lup[3][1] = 3; lup[3][2] = lup[3][3] = -1;
    lup[4][0] = 2; lup[4][1] = 3; lup[4][2] = lup[4][3] = -1;
    lup[5][0] = 0; lup[5][1] = 2; lup[5][2] = lup[5][3] = -1;
    lup[6][0] = 0; lup[6][1] = 1; lup[6][2] = 2; lup[6][3] = 3;
    lup[7][0] = 1; lup[7][1] = 2; lup[7][2] = lup[7][3] = -1;
    lup[8][0] = 1; lup[8][1] = 2; lup[8][2] = lup[8][3] = -1;
    lup[9][0] = 0; lup[9][1] = 3; lup[9][2] = 1; lup[9][3] = 2;
    lup[10][0] = 0; lup[10][1] = 2; lup[10][2] = lup[10][3] = -1;
    lup[11][0] = 2; lup[11][1] = 3; lup[11][2] = lup[11][3] = -1;
    lup[12][0] = 1; lup[12][1] = 3; lup[12][2] = lup[12][3] = -1;
    lup[13][0] = 0; lup[13][1] = 1; lup[13][2] = lup[13][3] = -1;
    lup[14][0] = 0; lup[14][1] = 3; lup[14][2] = lup[14][3] = -1;
    lup[15][0] = lup[15][1] = lup[15][2] = lup[15][3] = -1;

    // Evaluate isolines at each cell and call sl.AddSegment()
    float isovalue = 3.2;
    for(int i=0; i < dims[0] - 1; i++){
        for(int j=0; j < dims[1] - 1; j++){
            // Find ID's for vertex 3 of the current cell
            int logicalVertexId[2][2];
            int vertexId[2];
            logicalVertexId[0][0] = i + 1;
            logicalVertexId[0][1] = j + 1;
            vertexId[0] = GetPointIndex(logicalVertexId[0], dims);
            
            // Determine which case the cell fits
            int icase = 0;
            if(F[vertexId[0]] >= isovalue)
                icase += 8;
            logicalVertexId[0][0] -= 1;
            vertexId[0] = GetPointIndex(logicalVertexId[0], dims);
            if(F[vertexId[0]] >= isovalue)
                icase += 4;
            logicalVertexId[0][0] += 1;
            logicalVertexId[0][1] -= 1;
            vertexId[0] = GetPointIndex(logicalVertexId[0], dims);
            if(F[vertexId[0]] >= isovalue)
                icase += 2;
            logicalVertexId[0][0] -= 1;
            vertexId[0] = GetPointIndex(logicalVertexId[0], dims);
            if(F[vertexId[0]] >= isovalue)
                icase += 1;

            // Determine how many segments the cell will contain (0, 1, or 2)
            int nsegments = numSegments[icase];
            for(int k=0; k < nsegments; k++){
                float myX[2];
                float myY[2];

                // Find the endpoint locations for line segment. Runs twice for each segment
                for(int l=0; l < 2; l++){
                    // Find vertex ID's for the vertices to be interpolated 
                    bool horizontal;
                    switch(lup[icase][2*k+l]){
                        case 0: logicalVertexId[0][0] = i;
                            logicalVertexId[0][1] = j;
                            logicalVertexId[1][0] = i + 1;
                            logicalVertexId[1][1] = j;
                            horizontal = true;
                            break;
                        case 1: logicalVertexId[0][0] = i + 1;
                            logicalVertexId[0][1] = j;
                            logicalVertexId[1][0] = i + 1;
                            logicalVertexId[1][1] = j + 1;
                            horizontal = false;
                            break;
                        case 2: logicalVertexId[0][0] = i;
                            logicalVertexId[0][1] = j + 1;
                            logicalVertexId[1][0] = i + 1;
                            logicalVertexId[1][1] = j + 1;
                            horizontal = true;
                            break;
                        case 3: logicalVertexId[0][0] = i;
                            logicalVertexId[0][1] = j;
                            logicalVertexId[1][0] = i;
                            logicalVertexId[1][1] = j + 1;
                            horizontal = false;
                            break;
                        default: cerr << "ERROR!" << endl;
                    }
                    vertexId[0] = GetPointIndex(logicalVertexId[0], dims);
                    vertexId[1] = GetPointIndex(logicalVertexId[1], dims);

                    // Determine t-value for interpolation and solve the interpolation
                    float tValue = (isovalue - F[vertexId[0]]) / (F[vertexId[1]] - F[vertexId[0]]);
                    if(horizontal){
                        myX[l] = X[logicalVertexId[0][0]] + (tValue * (X[logicalVertexId[1][0]] - X[logicalVertexId[0][0]]));
                        myY[l] = Y[logicalVertexId[0][1]];
                    } else {
                        myX[l] = X[logicalVertexId[0][0]];
                        myY[l] = Y[logicalVertexId[0][1]] + (tValue * (Y[logicalVertexId[1][1]] - Y[logicalVertexId[0][1]]));
                    }
                }
                sl.AddSegment(myX[0], myY[0], myX[1], myY[1]);
            }
        }
    }

// MY CODE TO GENERATE ISOLINES ENDS HERE!


    vtkPolyData *pd = sl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
