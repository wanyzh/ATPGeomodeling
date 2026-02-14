// ATPExtrudeWithConstrain.cpp : This file contains the "main" function. Program execution begins and ends here.
//

#include <iostream>
#include "CATPGeoModel.h"
#include "vtkGAMBITReader.h"
#include "vtkSTLReader.h"
#include "vtkXMLMultiBlockDataWriter.h"
#include "vtkUnstructuredGridWriter.h"
double f(double a,double x, double y,double k,double phi)
{
	return  a * (sin(k*x + phi) + cos(k*x + phi));
}
int main()
{
	CATPGeoModel geomodel;
	
	/******************************************************************/
	//Read node indices for all fractures on the reference surface. Each fracture includes azimuth,
	//intersections with other fractures, and intersections with the boundary.
	//This surface is the reference for all other surfaces.
	vtkNew<vtkGAMBITReader> gambitReader;
	gambitReader->SetFileName("surf_frac_3.neu");
	gambitReader->Update();
	vtkUnstructuredGrid* refGrid = gambitReader->GetOutput();

	vtkSmartPointer<vtkUnstructuredGrid> surf1Grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	surf1Grid->DeepCopy(refGrid);

	auto points=refGrid->GetPoints();
	for (auto i=0;i<points->GetNumberOfPoints();++i)
	{
		double* coords = points->GetPoint(i);
		//coords[2] -= 0.5; // Update Z coordinate
		points->SetPoint(i, coords);
	}
	
	vtkNew<vtkUnstructuredGridWriter> uGridWriter;
	uGridWriter->SetInputData(refGrid);
	uGridWriter->SetFileName("surf_frac_3.vtk");
	uGridWriter->Update();
	//return 0;	

	//Add small undulations to the surface Z to mimic real strata; real data does not need this step.
	
	auto updatedPoints = surf1Grid->GetPoints();
	for (int i = 0; i < surf1Grid->GetNumberOfPoints(); ++i)
	{
		double* coords = updatedPoints->GetPoint(i);
		coords[2] += (sin(coords[0] * 0.1) + cos(coords[1] * 0.1)); // Update Z coordinate
		updatedPoints->SetPoint(i, coords);
	}
	//vtkNew<vtkUnstructuredGridWriter> uGridWriter;
	uGridWriter->SetInputData(surf1Grid);
	uGridWriter->SetFileName("surf1Grid.vtk");
	uGridWriter->Update();


	//Set polygon node indices in surf1Grid (sorted pseudo-clockwise)
	//bound boundPI;
	//boundPI._p.emplace_back(49); boundPI._p.emplace_back(69); boundPI._p.emplace_back(59); boundPI._p.emplace_back(94);


	//Populate node indices, azimuth, and intersections based on input fractures.
	//This simulates how the main program fills fracture data.
	
	int numFrac = 3;//3 fractures
	std::vector<frac> vecFracs;
	vecFracs.resize(numFrac);
	
	//Fracture 1 (two intersections with the boundary)
	vecFracs[0].p1 = 127; vecFracs[0].p2 = 109; vecFracs[0].beta = M_PI / 12;
	//Whether the fracture endpoints lie on the boundary
	vecFracs[0].isOnBound[0] = 0; vecFracs[0].isOnBound[1] = 0;
	//Intersections with other fractures (intersetionPoint): the first value is the node index in surf1Grid,
	//the second value is the other fracture's index.
	vecFracs[0].vecFracIntersect.emplace_back(intersetionPoint(1, 1)); vecFracs[0].vecFracIntersect.emplace_back(intersetionPoint(0, 2));
	//Fracture mesh node indices
	vector<int> temp = { 127, 133, 132, 131, 130, 129, 128, 1, 8, 7, 6, 5, 4, 3, 2, 0, 114, 113, 112, 111, 110, 109 };
	vecFracs[0].vecNodes.insert(vecFracs[0].vecNodes.end(),temp.begin(),temp.end());

	//Fracture 2
	vecFracs[1].p1 = 115; vecFracs[1].p2 = 121; vecFracs[1].beta = M_PI / 12;
	vecFracs[1].vecFracIntersect.emplace_back(intersetionPoint(1, 0)); 
	vector<int> temp2 = { 115,120,119,118,117,116,1,126,125,124,123,122,121 };
	vecFracs[1].vecNodes.insert(vecFracs[1].vecNodes.end(), temp2.begin(),temp2.end());
	
	//Fracture 3
	vecFracs[2].p1 = 134; vecFracs[2].p2 = 139; vecFracs[2].beta = M_PI / 12;
	vecFracs[2].vecFracIntersect.emplace_back(intersetionPoint(0, 0));
	vector<int> temp3 = { 134,135,136,137,138,0,140,141,142,139};
	vecFracs[2].vecNodes.insert(vecFracs[2].vecNodes.end(), temp3.begin(),temp3.end());
	


	
	/******************************************************************/
	//Store all surfaces in a vtkMultiBlockDataSet
	vtkSmartPointer<vtkMultiBlockDataSet> allOriSurfsGrid = vtkSmartPointer<vtkMultiBlockDataSet>::New();
	allOriSurfsGrid->SetBlock(0, surf1Grid);

	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName("surf2.stl");
	reader->Update();
	vtkSmartPointer<vtkUnstructuredGrid> surf2Grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	surf2Grid->DeepCopy(reader->GetOutput());
	allOriSurfsGrid->SetBlock(1, surf2Grid);

	vtkSmartPointer<vtkSTLReader> reader1 = vtkSmartPointer<vtkSTLReader>::New();
	reader1->SetFileName("surf3.stl");
	reader1->Update();
	vtkSmartPointer<vtkUnstructuredGrid> surf3Grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	surf3Grid->DeepCopy(reader1->GetOutput());
	allOriSurfsGrid->SetBlock(2, surf3Grid);
	

	vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
	writer->SetInputData(allOriSurfsGrid);
	writer->SetFileName("allOriSurfsGrid.vtm");
	writer->Write();

	/******************************************************************/
	//Pass fractures and surfaces to CATPGeoModel for processing
	geomodel.SetInputData(allOriSurfsGrid, vecFracs, refGrid);
	//geomodel.SetInputData(allOriSurfsGrid, refGrid);
	//geomodel.ImportVecFracs("fracs.txt");

	auto grid = geomodel.Generate();
	geomodel.Export2VtkFile("geomodelGrid.vtk", grid);
	geomodel.Export2VtmFile("surfOutput.vtm", geomodel.GetSurfOutPut());
	geomodel.Export2VtmFile("geomodelFracs.vtm", geomodel.GetFracs());

	vector<int>  numOfSubLayers{ 4, 10 };
	geomodel.Export2VtkFile("geomodelFineGrid.vtk", geomodel.RefineLayer(numOfSubLayers));
	geomodel.Export2VtmFile("geomodelFineFracs.vtm", geomodel.GetFracs());
	geomodel.Export2VtmFile("surfRefined.vtm", geomodel.getSurfRefineOutputMB());
	return 0;
}

