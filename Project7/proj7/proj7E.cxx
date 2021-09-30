
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

   vtkDataSetMapper *mapper1 = vtkDataSetMapper::New();
   mapper1->SetInputConnection(cut->GetOutputPort());

   vtkDataSetMapper *mapper2 = vtkDataSetMapper::New();
   mapper2->SetInputConnection(cf->GetOutputPort());
   
   vtkLookupTable *lut = vtkLookupTable::New();
   mapper1->SetLookupTable(lut);
   mapper1->SetScalarRange(1,6);
   mapper2->SetLookupTable(lut);
   mapper2->SetScalarRange(1,6);
   for(int i = 0; i < 256; i++){
      lut->SetTableValue(i, i, 0, (255 - i));
   }
   lut->Build();

   vtkActor *actor1 = vtkActor::New();
   actor1->SetMapper(mapper1);

   vtkActor *actor2 = vtkActor::New();
   actor2->SetMapper(mapper2);

   vtkRenderer *ren = vtkRenderer::New();
   ren->AddActor(actor1);
   ren->AddActor(actor2);

   vtkRenderWindow *renwin = vtkRenderWindow::New();
   renwin->AddRenderer(ren);
   renwin->SetSize(768, 768);

   vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
   iren->SetRenderWindow(renwin);
   renwin->Render();
   iren->Start();
}


