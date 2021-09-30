
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkContourFilter.h>

int main(int argc, char *argv[])
{
   vtkDataSetReader *reader = vtkDataSetReader::New();
   reader->SetFileName("noise.vtk");
   reader->Update();

   vtkPlane *plane = vtkPlane::New();
   plane->SetNormal(0,0,1);

   vtkCutter *cut = vtkCutter::New();
   cut->SetCutFunction(plane);
   cut->SetInputConnection(reader->GetOutputPort());

   vtkContourFilter *cf = vtkContourFilter::New();
   cf->SetNumberOfContours(2);
   cf->SetValue(0, 2.4);
   cf->SetValue(1, 4.0);
   cf->SetInputConnection(reader->GetOutputPort());

   vtkDataSetMapper *mapperLeft = vtkDataSetMapper::New();
   mapperLeft->SetInputConnection(cut->GetOutputPort());

   vtkDataSetMapper *mapperRight = vtkDataSetMapper::New();
   mapperRight->SetInputConnection(cf->GetOutputPort());
   
   vtkLookupTable *lut = vtkLookupTable::New();
   mapperLeft->SetLookupTable(lut);
   mapperLeft->SetScalarRange(1,6);
   mapperRight->SetLookupTable(lut);
   mapperRight->SetScalarRange(1,6);
   for(int i = 0; i < 256; i++){
      lut->SetTableValue(i, i, 0, (255 - i));
   }
   lut->Build();

   vtkActor *actorLeft = vtkActor::New();
   actorLeft->SetMapper(mapperLeft);

   vtkActor *actorRight = vtkActor::New();
   actorRight->SetMapper(mapperRight);

   vtkRenderer *renLeft = vtkRenderer::New();
   renLeft->AddActor(actorLeft);
   renLeft->SetViewport(0.0, 0.0, 0.5, 1.0);

   vtkRenderer *renRight = vtkRenderer::New();
   renRight->AddActor(actorRight);
   renRight->SetViewport(0.5, 0.0, 1.0, 1.0);

   vtkRenderWindow *renwin = vtkRenderWindow::New();
   renwin->AddRenderer(renLeft);
   renwin->AddRenderer(renRight);
   renwin->SetSize(768, 768);

   vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
   iren->SetRenderWindow(renwin);
   renwin->Render();
   iren->Start();
}
