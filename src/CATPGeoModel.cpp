#include "CATPGeoModel.h"
#include <vtkUnstructuredGridReader.h>
#include <vtkGeometryFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointData.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vtkPointLocator.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXmlUnstructuredGridWriter.h>
#include <vtkGAMBITReader.h>
#include <vtkCellLocator.h>
#include <vtkDataSet.h>
#include <vtkSTLReader.h>
#include <vtkLine.h>
#include "vtkXMLMultiBlockDataWriter.h"
#include <vtkCelldata.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPlane.h>
//#include <Eigen/PardisoSupport>
#include <vtkPlaneSource.h>
#include "vtkCutter.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyData.h"
#include "vtkXMLPolyDataWriter.h"
#include <algorithm> // Includes std::max and std::min
#include <iostream>
typedef Eigen::SparseMatrix<double> SpMat; // Sparse matrix type
typedef Eigen::VectorXd Rsv; // Vector type

using namespace std;
//#define PI 3.1415926;
CATPGeoModel::CATPGeoModel()
{
	m_fracsGrid = vtkSmartPointer<vtkMultiBlockDataSet>::New();
	m_surfOutputMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
	m_surfRefineOutputMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
	m_dMaxLength = 500;
}
void CATPGeoModel::SetInputData(vtkMultiBlockDataSet* allSurfs, vector<frac> fracs,vtkUnstructuredGrid* refGrid)
{	
	m_surfInputMB = allSurfs;
	m_vecFracs = fracs;
	m_referenceGrid = refGrid;
}
void CATPGeoModel::SetInputData(vtkMultiBlockDataSet* allSurfs,  vtkUnstructuredGrid* refGrid)
{
	m_surfInputMB = allSurfs;	
	m_referenceGrid = refGrid;
}
vtkSmartPointer<vtkUnstructuredGrid> CATPGeoModel::Generate()
{	
	ExportVecFracs("fracs.txt");
	/************************************************************************/
	//Build rays from fracture endpoints (longitudinal boundary lines). No fracture nodes overlap.
	vector<vector<frac>> rayFracs;//Two rays emitted from each fracture's endpoints
	for (int i = 0; i < m_vecFracs.size(); ++i)
	{	
		//Update vecPoints with node coordinates on the fracture
		m_vecFracs[i].vecPoints.resize(m_vecFracs[i].vecNodes.size());
		for (int j = 0; j < m_vecFracs[i].vecNodes.size(); ++j)
		{
			m_referenceGrid->GetPoint(m_vecFracs[i].vecNodes[j], m_vecFracs[i].vecPoints[j]._x);
		}

		double alpha = m_vecFracs[i].beta;
		double p1[3], p2[3];
		m_referenceGrid->GetPoint(m_vecFracs[i].p1, p1);
		m_referenceGrid->GetPoint(m_vecFracs[i].p2, p2);

		// Define a point and direction vector of the original line
		//Eigen::Vector3d P0(p1[0], p1[1], p1[2]);
		Eigen::Vector3d V(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
		//Eigen::Vector3d V1(p2[0] - p1[0], p2[1] - p1[1] + 10, p2[2] - p1[2]);
		V.normalize(); // Normalize direction vector
		//V1.normalize();

		// Compute a direction vector perpendicular to the original line
		//Eigen::Vector3d direction_perpendicular = V.cross(V1);
		Eigen::Vector3d direction_perpendicular(0, 0, 1);//Reference surface normal is the Z axis
		direction_perpendicular.normalize();

		//Rotate direction vector
		Eigen::AngleAxisd rotation(alpha,V);//Rotate around V
		Eigen::Vector3d rotate_direction = rotation * direction_perpendicular;
		//double dLength = m_dMaxLength;//Segment length
		//Rays from fracture endpoints
		frac f1, f2;

		f1.index = m_vecFracs[i].p1;
		f1.frac_index = i;
		f2.index = m_vecFracs[i].p2;
		f2.frac_index = i;

		f1.coP1[0] = p1[0] - m_dMaxLength * rotate_direction(0);
		f1.coP1[1] = p1[1]- m_dMaxLength * rotate_direction(1);
		f1.coP1[2] = p1[2]- m_dMaxLength * rotate_direction(2);
		
		f1.coP2[0] = p1[0] + m_dMaxLength * rotate_direction(0);
		f1.coP2[1] = p1[1] + m_dMaxLength * rotate_direction(1);
		f1.coP2[2] = p1[2] + m_dMaxLength * rotate_direction(2);

		f2.coP1[0] = p2[0]- m_dMaxLength * rotate_direction(0);
		f2.coP1[1] = p2[1]- m_dMaxLength * rotate_direction(1);
		f2.coP1[2] = p2[2]- m_dMaxLength * rotate_direction(2);
		
		f2.coP2[0] = p2[0] + m_dMaxLength * rotate_direction(0);
		f2.coP2[1] = p2[1] + m_dMaxLength * rotate_direction(1);
		f2.coP2[2] = p2[2] + m_dMaxLength * rotate_direction(2);

		vector<frac> tempFrac;		
		tempFrac.emplace_back(f1);		
		tempFrac.emplace_back(f2);
		rayFracs.emplace_back(tempFrac);		
	}

	/************************************************************************/
	//Build rays from intersections between fractures; these rays must account for both dip angles.
	//Compute the intersection line by intersecting two planes.
	//For two fracture planes: compute their normals, and the cross product gives the intersection direction.
	//For direction consistency, enforce all directions to point towards +Z.
	vector<vector<frac>> rayFracsInterWithFracs;//Two rays emitted from each fracture-fracture intersection
	for (int i=0;i<m_vecFracs.size();++i)
	{
		vector<frac> tempFrac;
		for (int j=0;j<m_vecFracs[i].vecFracIntersect.size();++j)
		{
			double p1[3], p2[3];
			m_referenceGrid->GetPoint(m_vecFracs[i].p1, p1);
			m_referenceGrid->GetPoint(m_vecFracs[i].p2, p2);

			Eigen::Vector3d V(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
			Eigen::Vector3d V1(rayFracs[i][0].coP2[0] - p1[0], rayFracs[i][0].coP2[1] - p1[1] , rayFracs[i][0].coP2[2] - p1[2]);
			Eigen::Vector3d normal1 = V.cross(V1);						

			int fracIndex = m_vecFracs[i].vecFracIntersect[j].index;
			double pj1[3], pj2[3];
			m_referenceGrid->GetPoint(m_vecFracs[fracIndex].p1, pj1);
			m_referenceGrid->GetPoint(m_vecFracs[fracIndex].p2, pj2);

			Eigen::Vector3d Vj(pj2[0] - pj1[0], pj2[1] - pj1[1], pj2[2] - pj1[2]);
			Eigen::Vector3d Vj1(rayFracs[fracIndex][0].coP2[0] - pj1[0], rayFracs[fracIndex][0].coP2[1] - pj1[1], rayFracs[fracIndex][0].coP2[2] - pj1[2]);
			Eigen::Vector3d normal2 = Vj.cross(Vj1);

			Eigen::Vector3d crossNormal = normal1.cross(normal2);//Intersection line direction vector
			crossNormal.normalize();
			double angle = std::acos(crossNormal.z()); // Angle to +Z

			if (angle > M_PI / 2) { // If the angle is greater than 90 degrees
				crossNormal = -crossNormal; // Flip direction
			}

			double dLength = m_dMaxLength;//Segment length

			frac f1;
			f1.index = m_vecFracs[i].vecFracIntersect[j].pI;
			f1.frac_index = i;

			double tempP[3];
			m_referenceGrid->GetPoint(m_vecFracs[i].vecFracIntersect[j].pI, tempP);
			
			f1.coP1[0]= tempP[0] - dLength * crossNormal(0);
			f1.coP1[1] = tempP[1] - dLength * crossNormal(1);
			f1.coP1[2] = tempP[2] - dLength * crossNormal(2);
			
			f1.coP2[0] = tempP[0] + dLength * crossNormal(0);
			f1.coP2[1] = tempP[1] + dLength * crossNormal(1);
			f1.coP2[2] = tempP[2] + dLength * crossNormal(2);
			tempFrac.emplace_back(f1);
		}
		rayFracsInterWithFracs.emplace_back(tempFrac);
	}

	lineGridConsructAndOutput(rayFracs, "fraFracs.vtk");
	lineGridConsructAndOutput(rayFracsInterWithFracs, "rayFracsInterWithFracs.vtk");
	m_nPointsNumberOfFracs = rayFracs.size();	

	/************************************************************************/
	vector<vector<frac>> allNewFracs;//All fracture data on all surfaces (endpoints, intersection points, etc.)
	for (int j = 0; j < m_surfInputMB->GetNumberOfBlocks(); ++j)
	{
		vector<frac> fracs;//All fractures on surface j
		for (int i = 0; i < m_vecFracs.size(); ++i)
		{			
			frac frac_j(m_vecFracs[i]);//Fracture i on surface j
			frac_j.isOnBound[0] = 0; frac_j.isOnBound[1] = 0;//Initialize: endpoints not on boundary (updated later)
			/**************************************************/
			//Compute fracture endpoints on the new surface.
			//Endpoints are obtained by intersecting endpoint rays with the new surface. Cases:
			//1) Intersection lies inside the new surface: use it directly as the endpoint.
			//2) Intersection lies on the boundary or outside: compute the intersection line between
			//   the fracture plane and the new surface, and take its boundary points as endpoints.

			//Intersect endpoint rays with the new surface
			intersetionPoint tempP1,tempP2;//Intersection points
			vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
			auto grid = vtkUnstructuredGrid::SafeDownCast(m_surfInputMB->GetBlock(j));
			cellLocator->SetDataSet(grid);
			cellLocator->BuildLocator();
			double t;//Parametric coordinate of the closest point on the line
			double pcoords[3];//Parametric coordinates inside the cell (local coordinate system)
			vtkIdType cellId;//Intersected cell id
			int subId;//Intersected sub-cell id (if applicable)
			auto res1 = cellLocator->IntersectWithLine(rayFracs[i][0].coP1, rayFracs[i][0].coP2, 0.0, t, tempP1.p, pcoords, subId, cellId);
			auto res2 = cellLocator->IntersectWithLine(rayFracs[i][1].coP1, rayFracs[i][1].coP2, 0.0, t, tempP2.p, pcoords, subId, cellId);
			
			frac_j.vecPoints[0] = Point(tempP1.p);			
			frac_j.vecPoints[frac_j.vecPoints.size() - 1] = Point(tempP2.p);
			for (int tempi=0;tempi<3;++tempi)
			{
				frac_j.coP1[tempi] = tempP1.p[tempi];
				frac_j.coP2[tempi] = tempP2.p[tempi];
			}

			//If the intersection is not inside the new surface, compute the fracture-plane / surface
			//intersection curve and use its boundary points as fracture endpoints on this surface.
			if (0==res1 || 0==res2)
			{
				double p1[3], p2[3];
				m_referenceGrid->GetPoint(m_vecFracs[i].p1, p1);
				m_referenceGrid->GetPoint(m_vecFracs[i].p2, p2);
				Eigen::Vector3d V(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
				Eigen::Vector3d V1(rayFracs[i][0].coP2[0] - p1[0], rayFracs[i][0].coP2[1] - p1[1], rayFracs[i][0].coP2[2] - p1[2]);
				Eigen::Vector3d normal1 = V.cross(V1);
				normal1.normalize();

				// Create a cutting plane
				vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
				plane->SetOrigin(p1); // Set plane origin
				plane->SetNormal(normal1(0), normal1(1), normal1(2)); // Set plane normal
				//The fracture plane is defined by its normal


				// Create a cutter and set the cutting plane
				vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
				auto grid = vtkUnstructuredGrid::SafeDownCast(m_surfInputMB->GetBlock(j));
				cutter->SetInputData(grid);
				cutter->SetCutFunction(plane);

				// Write the result to file
				// Generate file name
				std::stringstream ss;
				ss << "output_" << j << ".vtp"; // Name based on loop index
				vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
				writer->SetFileName(ss.str().c_str()); // Set output file name
				writer->SetInputConnection(cutter->GetOutputPort());
				writer->Write();

				//Get start and end coordinates of the cut curve.
				//The curve points are not ordered, so sort points to determine endpoints.
				std::vector<Point> tempVecPoints;
				int tempNum = cutter->GetOutput()->GetNumberOfPoints();
				for (int ci = 0; ci < tempNum; ++ci)
				{
					tempVecPoints.emplace_back(Point(cutter->GetOutput()->GetPoint(ci)));
				}
				sortPointsByPolarAngle(tempVecPoints, Point(0, 0, -1000));

				//Output sorted points
				//vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
				//vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
				//for (int vp=0;vp< cutter->GetOutput()->GetNumberOfPoints();++vp)
				//{
				//	points->InsertNextPoint(tempVecPoints[vp]._x);
				//}
				//polyData->SetPoints(points);				
				//vtkSmartPointer <vtkXMLPolyDataWriter> polywriter = vtkSmartPointer <vtkXMLPolyDataWriter>::New();
				//polywriter->SetFileName("sorted_points.vtp"); // Set output file name
				//polywriter->SetInputData(polyData);
				//polywriter->Write();

				//Assign sorted endpoints back to the fracture.
				//To decide which endpoint corresponds to the original fracture: connect endpoints to form two lines,
				//and test whether the line intersection lies inside the quadrilateral formed by the four points.

				//Compute line intersection
				auto intersection = calculateIntersection(tempVecPoints[0], Point(p1), tempVecPoints[tempNum - 1], Point(p2));

				//If the two lines have no intersection, or the intersection lies outside the quadrilateral
				//formed by the four endpoints, the chosen endpoint pairing is considered correct.
				if (std::isnan(intersection._x[0]) || !isInside(tempVecPoints[0], Point(p1), tempVecPoints[tempNum - 1], Point(p2), intersection))
				{
					if (0==res1 || (res1 && m_vecFracs[i].isOnBound[0])) //If the intersection is not inside the surface, or the reference endpoint is on boundary, clamp endpoint to boundary
					{
						frac_j.isOnBound[0] = 1;
						frac_j.vecPoints[0] = tempVecPoints[0];
						for (int tempi = 0; tempi < 3; ++tempi)
						{
							frac_j.coP1[tempi] = tempVecPoints[0]._x[tempi];
						}
					}					
					if (0==res2 || (res2 && m_vecFracs[i].isOnBound[1])) {
						frac_j.isOnBound[1] = 1;
						frac_j.vecPoints[frac_j.vecPoints.size() - 1] = tempVecPoints[tempNum - 1];
						for (int tempi = 0; tempi < 3; ++tempi)
						{
							frac_j.coP2[tempi] = tempVecPoints[tempNum - 1]._x[tempi];
						}
					}					
				}
				else
				{
					if (0==res1 || (res1 && m_vecFracs[i].isOnBound[0])) {
						frac_j.isOnBound[0] = 1;
						frac_j.vecPoints[0] = tempVecPoints[tempNum - 1];
						for (int tempi = 0; tempi < 3; ++tempi)
						{
							frac_j.coP1[tempi] = tempVecPoints[tempNum - 1]._x[tempi];
						}
					}					
					if (0==res2 || ( res2 &&m_vecFracs[i].isOnBound[1] )) {
						frac_j.isOnBound[1] = 1;
						frac_j.vecPoints[frac_j.vecPoints.size() - 1] = tempVecPoints[0];
						for (int tempi = 0; tempi < 3; ++tempi)
						{
							frac_j.coP2[tempi] = tempVecPoints[0]._x[tempi];
						}
					}					
				}
			
			}

			/**************************************************/
			//Compute fracture-fracture intersection points on the new surface
			//The intersection between the inter-fracture line and the new surface may lie outside; ignore for now.
			frac_j.vecFracIntersect.resize(m_vecFracs[i].vecFracIntersect.size());
			for(int it=0;it<rayFracsInterWithFracs[i].size();++it)
			{
				intersetionPoint tempPi(m_vecFracs[i].vecFracIntersect[it]);				
				auto resi = cellLocator->IntersectWithLine(rayFracsInterWithFracs[i][it].coP1, rayFracsInterWithFracs[i][it].coP2, 0.0, t, tempPi.p, pcoords, subId, cellId);
				frac_j.vecFracIntersect[it] = tempPi;
				int nodeIndex = std::find(frac_j.vecNodes.begin(), frac_j.vecNodes.end(), tempPi.pI) - frac_j.vecNodes.begin();//Index of the intersection within the fracture node list
				frac_j.vecPoints[nodeIndex] = Point(tempPi.p);
			}

			/**************************************************/
			//Compute coordinates of the other fracture nodes on the new surface
			//Redistribute positions based on distances between endpoints and intersections
			//After building rays, compute coordinates via intersection
			std::vector<int> tempKeyIndex;
			tempKeyIndex.emplace_back(0);			
			for (int oi=0;oi< m_vecFracs[i].vecFracIntersect.size();++oi)
			{
				int nodeIndex = std::find(frac_j.vecNodes.begin(), frac_j.vecNodes.end(), m_vecFracs[i].vecFracIntersect[oi].pI) - frac_j.vecNodes.begin();
				tempKeyIndex.emplace_back(nodeIndex);
			}
			tempKeyIndex.emplace_back(frac_j.vecNodes.size()-1);
			for (int it=0;it<tempKeyIndex.size()-1;++it)
			{
				double length = distanceTwoPoints(m_vecFracs[i].vecPoints[tempKeyIndex[it]], m_vecFracs[i].vecPoints[tempKeyIndex[it + 1]]);
				double lengthBar = distanceTwoPoints(frac_j.vecPoints[tempKeyIndex[it]], frac_j.vecPoints[tempKeyIndex[it + 1]]);
				for (int jt= tempKeyIndex[it];jt<tempKeyIndex[it+1];++jt)
				{
					double Lj = distanceTwoPoints(m_vecFracs[i].vecPoints[tempKeyIndex[it]], m_vecFracs[i].vecPoints[jt]);
					double ratio = Lj / length;
					
					double tempPjt[3]{
						frac_j.vecPoints[tempKeyIndex[it]]._x[0] + ratio * (frac_j.vecPoints[tempKeyIndex[it + 1]]._x[0] - frac_j.vecPoints[tempKeyIndex[it]]._x[0]) ,
						frac_j.vecPoints[tempKeyIndex[it]]._x[1] + ratio * (frac_j.vecPoints[tempKeyIndex[it + 1]]._x[1] - frac_j.vecPoints[tempKeyIndex[it]]._x[1]),
						frac_j.vecPoints[tempKeyIndex[it]]._x[2] + ratio * (frac_j.vecPoints[tempKeyIndex[it + 1]]._x[2] - frac_j.vecPoints[tempKeyIndex[it]]._x[2]) };

					Eigen::Vector3d V(tempPjt[0] - m_vecFracs[i].vecPoints[jt]._x[0], tempPjt[1] - m_vecFracs[i].vecPoints[jt]._x[1], tempPjt[2] - m_vecFracs[i].vecPoints[jt]._x[2]);

					double p1[3];
					double p2[3];

					p1[0] = tempPjt[0]-m_dMaxLength * V(0);
					p1[1] = tempPjt[1] -m_dMaxLength * V(1);
					p1[2] = tempPjt[2] -m_dMaxLength * V(2);

					p2[0] = /*m_vecFracs[i].vecPoints[jt]._x[0];//*/tempPjt[0] +m_dMaxLength * V(0);
					p2[1] = /*m_vecFracs[i].vecPoints[jt]._x[1];//*/tempPjt[1] +m_dMaxLength * V(1);
					p2[2] = /*m_vecFracs[i].vecPoints[jt]._x[2];//*/tempPjt[2] +m_dMaxLength * V(2);

					double tempP[3];
					auto resi = cellLocator->IntersectWithLine(p1, p2, 1e-5, t, tempP, pcoords, subId, cellId);

					frac_j.vecPoints[jt]._x[0] = tempP[0];
					frac_j.vecPoints[jt]._x[1] = tempP[1];
					frac_j.vecPoints[jt]._x[2] = tempP[2];

				}				
			}
			fracs.emplace_back(frac_j);
			
			
		}
		allNewFracs.emplace_back(fracs);	
		/*************************************/
		//Conformal parameterization to obtain a planar mesh
		//First collect all fixed points on the new surface and record their reference indices
		//Fixed points include all points on fractures
		vector<intersetionPoint> keyPoints;
		for (int i = 0; i < fracs.size(); ++i)
		{
			for (int k=0;k<fracs[i].vecNodes.size();++k)
			{
				keyPoints.emplace_back(intersetionPoint(fracs[i].vecPoints[k],fracs[i].vecNodes[k]));
			}
		}
		//Conformal parameterization
		auto grid = Floater(m_referenceGrid, keyPoints);//keyPoints are the fracture constraints for parameterization
		m_surfOutputMB->SetBlock(j, grid);
	}
	//fracGridConsructAndOutput(allNewFracs, "allNewFracs.vtm");
	//Export2VtmFile("surfOutputMBBefore.vtm", m_surfOutputMB);
	/*************************************/
	//Interpolate the parameterized surfaces back to match real surface undulations.
	//Project each parameterized point along the vertical direction onto the real surface and take its Z.
	for (int i = 0; i < m_surfOutputMB->GetNumberOfBlocks(); ++i)
	{
		// Create a vtkCellLocator and set vtkUnstructuredGrid as the source
		auto grid1 = vtkUnstructuredGrid::SafeDownCast(m_surfInputMB->GetBlock(i));
		vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
		cellLocator->SetDataSet(grid1);
		cellLocator->BuildLocator();

		auto grid = vtkUnstructuredGrid::SafeDownCast(m_surfOutputMB->GetBlock(i));
		for (int j = 0; j < grid->GetNumberOfPoints(); ++j)
		{

			// Output point coordinate buffer
			double projectedPoint[3];
			// Project each point
			//Find the closest points to TestPoint
			auto assistCell = vtkSmartPointer<vtkGenericCell>::New();
			double closestPoint[3];//the coordinates of the closest point will be returned here
			double closestPointDist2; //the squared distance to the closest point will be returned here
			vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
			int subId;
			cellLocator->FindClosestPoint(grid->GetPoint(j), closestPoint, assistCell, cellId, subId, closestPointDist2);

			double* coords = grid->GetPoint(j);
			coords[2] = closestPoint[2];
			grid->GetPoints()->SetPoint(j, coords);
		}
	}

	/*************************************/
	//Connect all surface meshes by topology to build the 3D mesh

	vtkSmartPointer<vtkUnstructuredGrid> ATPGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> ATPPoints = vtkSmartPointer<vtkPoints>::New();

	int numPointsOfOneLayer = m_referenceGrid->GetNumberOfPoints();	
	for (int i = 0; i < m_surfOutputMB->GetNumberOfBlocks(); ++i)
	{
		auto grid = vtkUnstructuredGrid::SafeDownCast(m_surfOutputMB->GetBlock(i));
		for (int j = 0; j < grid->GetNumberOfPoints(); ++j)
		{
			ATPPoints->InsertNextPoint(grid->GetPoint(j));
		}
	}
	ATPGrid->SetPoints(ATPPoints);

	//Build cells
	//Create and set cell scalar data: the layer index for each cell
	vtkSmartPointer<vtkIdTypeArray> cellLayer = vtkSmartPointer<vtkIdTypeArray>::New();
	cellLayer->SetName("Layer");
	for (int i = 0; i < m_surfOutputMB->GetNumberOfBlocks()-1; ++i)
	{
		auto grid = vtkUnstructuredGrid::SafeDownCast(m_surfOutputMB->GetBlock(i));
		for (int j = 0; j < grid->GetNumberOfCells(); ++j)
		{
			auto cell = grid->GetCell(j);
			auto pIds = cell->GetPointIds();
			vtkIdType prism[6] = {
				pIds->GetId(0) + i * numPointsOfOneLayer, pIds->GetId(1) + i * numPointsOfOneLayer ,pIds->GetId(2) + i * numPointsOfOneLayer ,
				pIds->GetId(0) + (i + 1)*numPointsOfOneLayer,pIds->GetId(1) + (i + 1)*numPointsOfOneLayer ,pIds->GetId(2) + (i + 1)*numPointsOfOneLayer };

			ATPGrid->InsertNextCell(VTK_WEDGE, 6, prism);
			cellLayer->InsertNextValue(i);
		}
	}
	ATPGrid->GetCellData()->SetScalars(cellLayer);
	//Export2VtkFile("ATPGrid.vtk", ATPGrid);

	//Build fracture surface meshes
	fracsGridUpdate(m_surfOutputMB);
	
	//m_fracsGrid->SetPoints(FracPoints);
	return ATPGrid;
}
vtkSmartPointer<vtkMultiBlockDataSet> CATPGeoModel::GetFracs()
{
	return m_fracsGrid;
}
vtkSmartPointer<vtkMultiBlockDataSet> CATPGeoModel::GetSurfOutPut()
{
	return m_surfOutputMB;
}
vtkSmartPointer<vtkUnstructuredGrid> CATPGeoModel::RefineLayer(vector<int> numOfFineLayers)
{
	vtkSmartPointer<vtkUnstructuredGrid> ATPFineGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> ATPFinePoints = vtkSmartPointer<vtkPoints>::New();
	//vector<int> numOfFineLayers{ 4,5 };	//Insert intermediate layers between coarse layers (e.g., split each prism layer)
	vtkSmartPointer<vtkIdTypeArray> cellFineLayer = vtkSmartPointer<vtkIdTypeArray>::New();
	cellFineLayer->SetName("Layer");
	//auto surf1Grid = vtkUnstructuredGrid::SafeDownCast(m_surfOutputMB->GetBlock(0));
	
	int numPointsOfOneLayer = m_referenceGrid->GetNumberOfPoints();
	int tempIndex = -1;
	for (int i = 0; i < m_surfOutputMB->GetNumberOfBlocks() - 1; ++i)
	{	

		auto grid1 = vtkUnstructuredGrid::SafeDownCast(m_surfOutputMB->GetBlock(i));
		auto grid2 = vtkUnstructuredGrid::SafeDownCast(m_surfOutputMB->GetBlock(i + 1));
	
		for (int j = 0; j < grid1->GetNumberOfPoints(); ++j)
		{
			ATPFinePoints->InsertNextPoint(grid1->GetPoint(j));
		}
		tempIndex += 1;
		m_surfRefineOutputMB->SetBlock(tempIndex, grid1);

		for (int k = 1; k < numOfFineLayers[i]; ++k)//Iterate intermediate layers between grid1 and grid2
		{
			vtkSmartPointer<vtkUnstructuredGrid> refineLayerGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
			refineLayerGrid->DeepCopy(m_referenceGrid);
			auto rlPoints=refineLayerGrid->GetPoints();
			for (int l = 0; l < numPointsOfOneLayer; ++l)
			{
				double p1[3], p2[3];
				grid1->GetPoint(l, p1);
				grid2->GetPoint(l, p2);

				double totalLength = std::sqrt(std::pow(p2[0] - p1[0], 2) 
					+ std::pow(p2[1] - p1[1], 2) + std::pow(p2[2] - p1[2], 2));
				double stepSize = totalLength / (numOfFineLayers[i]);
				double t = k * stepSize / totalLength;
				double dividedPoint[3];
				dividedPoint[0] = p1[0] + (p2[0] - p1[0]) * t;
				dividedPoint[1] = p1[1] + (p2[1] - p1[1]) * t;
				dividedPoint[2] = p1[2] + (p2[2] - p1[2]) * t;
				rlPoints->SetPoint(l, dividedPoint);
				ATPFinePoints->InsertNextPoint(dividedPoint);
			}
			tempIndex += 1;
			m_surfRefineOutputMB->SetBlock(tempIndex , refineLayerGrid);
		}
		//for (int j = 0; j < grid2->GetNumberOfPoints(); ++j)
		//{
		//	ATPFinePoints->InsertNextPoint(grid2->GetPoint(j));
		//}
		//tempIndex += 1;
		//m_surfRefineOutputMB->SetBlock(tempIndex, grid2);

	}
	auto grid2 = vtkUnstructuredGrid::SafeDownCast(m_surfOutputMB->GetBlock(m_surfOutputMB->GetNumberOfBlocks() - 1));
	for (int j = 0; j < grid2->GetNumberOfPoints(); ++j)
	{
		ATPFinePoints->InsertNextPoint(grid2->GetPoint(j));
	}
	tempIndex += 1;
	m_surfRefineOutputMB->SetBlock(tempIndex, grid2);
	ATPFineGrid->SetPoints(ATPFinePoints);
	
	for (int i = 0; i < tempIndex; ++i)
	{
		for (int j = 0; j < m_referenceGrid->GetNumberOfCells(); ++j)
		{
			auto cell = m_referenceGrid->GetCell(j);
			auto pIds = cell->GetPointIds();
			vtkIdType prism[6] = {
				pIds->GetId(0) + i * numPointsOfOneLayer, pIds->GetId(1) + i * numPointsOfOneLayer ,pIds->GetId(2) + i * numPointsOfOneLayer ,
				pIds->GetId(0) + (i + 1)*numPointsOfOneLayer,pIds->GetId(1) + (i + 1)*numPointsOfOneLayer ,pIds->GetId(2) + (i + 1)*numPointsOfOneLayer };

			ATPFineGrid->InsertNextCell(VTK_WEDGE, 6, prism);
			cellFineLayer->InsertNextValue(i);
		}
	}
	ATPFineGrid->GetCellData()->SetScalars(cellFineLayer);
	//Export2VtkFile("geomodelFineGrid.vtk", ATPFineGrid);

	//When refining the mesh, refine fractures as well
	fracsGridUpdate(m_surfRefineOutputMB);
	return ATPFineGrid;
}
vtkSmartPointer<vtkMultiBlockDataSet> CATPGeoModel::GetSurfRefineOutPut()
{
	return m_surfRefineOutputMB;
}
double CATPGeoModel::getAngle(Point center, Point p) {
	double dx = p._x[0] - center._x[0];
	double dy = p._x[1] - center._x[1];
	return atan2(dy, dx);
}
bool  CATPGeoModel::comparePoints(Point center, Point a, Point b) {
	double angleA = getAngle(center, a);
	double angleB = getAngle(center, b);

	return angleA < angleB;
}
std::vector<Point> CATPGeoModel::sortPointsByPolarAngle(std::vector<Point>& points, Point center) {
	// Sort by polar angle around the center point
	std::sort(points.begin(), points.end(), [&](Point a, Point b) {
		return comparePoints(center, a, b);
	});

	return points;
}
double  CATPGeoModel::sumOfVectorElements(const std::vector<double>& vec) {
	double sum = 0.0;
	for (const auto& element : vec) {
		sum += element;
	}
	return sum;
}
bool CATPGeoModel::compareDoubleArrays(const double* array1, const double* array2, int size) {
	for (int i = 0; i < size; i++) {
		if (array1[i] != array2[i]) {
			return false;
		}
	}
	return true;
}
bool CATPGeoModel::compareTwoFrac(frac f1, frac f2)//Return true if the two fractures are equivalent
{
	if ((compareDoubleArrays(f1.coP1, f2.coP1, 3) || compareDoubleArrays(f1.coP1, f2.coP2, 3)) && (compareDoubleArrays(f1.coP2, f2.coP1, 3) || compareDoubleArrays(f1.coP2, f2.coP2, 3)))
	{
		return true;
	}
	else
		return false;
}
bool CATPGeoModel::IsfracInvectors(frac f, vector<frac> fracs)
{
	for (int i = 0; i < fracs.size(); ++i)
	{
		if (compareTwoFrac(f, fracs[i]))
		{
			return true;
		}
	}
	return false;
}
bool CATPGeoModel::IsPointAlreadyInserted(vtkPoints* points, double point[3])
{
	// Number of inserted points
	vtkIdType numPoints = points->GetNumberOfPoints();

	// Iterate inserted points
	for (vtkIdType i = 0; i < numPoints; i++)
	{
		double* existingPoint = points->GetPoint(i);
		// Check whether the point already exists
		if (fabs(existingPoint[0] - point[0]) < 1.0e-5 && fabs(existingPoint[1] - point[1]) < 1.0e-5 && fabs(existingPoint[2] - point[2]) < 1.0e-5)
		{
			return true;
		}
	}

	return false;
}
vtkUnstructuredGrid* CATPGeoModel::lineGridConsructAndOutput(vector<vector<frac>> vecFracs, string filename)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	for (int i = 0; i < vecFracs.size(); ++i)
	{
		for (int j=0;j< vecFracs[i].size();++j)
		{
			int index = i * (vecFracs[i].size() * 2) + j * 2;
			points->InsertPoint(index +0, vecFracs[i][j].coP1);
			points->InsertPoint(index +1, vecFracs[i][j].coP2);

			vtkSmartPointer<vtkLine> edge1 = vtkSmartPointer<vtkLine>::New();
			edge1->GetPointIds()->SetId(0, index);  // First point id
			edge1->GetPointIds()->SetId(1, index + 1);  // Second point id

			cells->InsertNextCell(edge1);
		}
		

	}
	grid->SetPoints(points);
	grid->SetCells(VTK_LINE, cells);

	vtkNew<vtkUnstructuredGridWriter> writer;
	writer->SetFileName(filename.c_str());
	writer->SetInputData(grid);
	writer->Update();

	return grid;

}

vtkMultiBlockDataSet* CATPGeoModel::fracGridConsructAndOutput(vector<vector<frac>> vecFracs, string filename)
{	
	vtkSmartPointer<vtkMultiBlockDataSet> allFracsGrid = vtkSmartPointer<vtkMultiBlockDataSet>::New();	
	int index = 0;
	for (int i = 0; i < vecFracs.size(); ++i)//Number of layers
	{
		for (int j = 0; j < vecFracs[i].size(); ++j)//Fractures per layer
		{
			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
			vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
			vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
			for (int k = 0; k < vecFracs[i][j].vecPoints.size()-1;++k)//Points per fracture
			{		

				points->InsertPoint(k, vecFracs[i][j].vecPoints[k]._x);	
				points->InsertPoint(k+1, vecFracs[i][j].vecPoints[k+1]._x);
				vtkSmartPointer<vtkLine> edge1 = vtkSmartPointer<vtkLine>::New();
				edge1->GetPointIds()->SetId(0, k);  // First point id
				edge1->GetPointIds()->SetId(1, k + 1);  // Second point id
				cells->InsertNextCell(edge1);
			}
			grid->SetPoints(points);
			grid->SetCells(VTK_LINE, cells);

			allFracsGrid->SetBlock(index, grid);
			++index;
		}
	}

	vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
	writer->SetInputData(allFracsGrid);
	writer->SetFileName(filename.c_str());
	writer->Write();
	return allFracsGrid;
	
}
std::vector<double> CATPGeoModel::Weigth(vtkPoints* points, vtkIdType node, std::vector<int> nodeAdjacents, WEIGHT_TYPE type)
{
	vtkNew<vtkPoints> vtk_points;
	vtk_points->InsertNextPoint(points->GetPoint(node));
	for (int i = 0; i < nodeAdjacents.size(); ++i)
	{
		vtk_points->InsertNextPoint(points->GetPoint(nodeAdjacents[i]));
	}
	// Create vtkPolyData and set point data
	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	polyData->SetPoints(vtk_points);

	
	int num = nodeAdjacents.size();
	std::vector<double> weight(num);
	switch (type)
	{
	case Tutte:
		for (int i = 0; i < num; ++i)
		{
			weight[i] = 1.0 / num;
		}
		break;
	case Floater_Shape_Preserve:
		break;
	case Floater_Mean_Value:
	{
		//Sort adjacent nodes counter-clockwise
		std::vector<Point> vec_points;
		for (int i = 0; i < num; ++i)
		{
			Point p(points->GetPoint(nodeAdjacents[i]), i);
			vec_points.emplace_back(p);
		}
		sortPointsByPolarAngle(vec_points, Point(points->GetPoint(node)));

		std::vector<double> vec_w;
		for (int i = 0; i < vec_points.size(); ++i)
		{
			int l1, l2;
			l1 = i - 1;
			l2 = i + 1;
			if (i == 0)
			{
				l1 = vec_points.size() - 1;

			}
			else if (i == vec_points.size() - 1)
			{
				l2 = 0;
			}
			Point PVi = vec_points[i] - Point(points->GetPoint(node));
			Point PVl1 = vec_points[l1] - Point(points->GetPoint(node));
			Point PVl2 = vec_points[l2] - Point(points->GetPoint(node));
			double alpha1 = acos(PVi * PVl1 / PVi.GetLength() / PVl1.GetLength());
			double alpha2 = acos(PVi * PVl2 / PVi.GetLength() / PVl2.GetLength());
			vec_w.emplace_back(tan(alpha1 / 2.0) + tan(alpha2 / 2.0) / PVi.GetLength());
		}
		double sum_w = sumOfVectorElements(vec_w);
		for (int i = 0; i < vec_points.size(); ++i)
		{
			int index = vec_points[i]._index;
			weight[index] = vec_w[i] / sum_w;
		}

	}
	break;
	case Discrete_Authalic_Parameterization:
		break;
	case Discrete_Conformal_Map:
		break;
	case Iterative_Authalic_Parameterization:
		break;
	default:
		break;
	}
	return weight;

}
int CATPGeoModel::findIndexInPoint(int index, vector<intersetionPoint> keyPoints)
{
	for (int i = 0; i < keyPoints.size(); ++i)
	{
		if (keyPoints[i].pI == index)
			return i;
	}
	return -1;
}

int CATPGeoModel::findIndexInPoint(int index, vector<Point> Points)
{
	for (int i = 0; i < Points.size(); ++i)
	{
		if (Points[i]._index == index)
			return i;
	}
	return -1;
}
int CATPGeoModel::findPoint(std::vector<Point> points, Point p)
{

	//Point p(reader->GetOutput()->GetPoint(i));
	auto it = std::find(points.begin(), points.end(), p); // Iterator to element in vector
	int index = -1;
	if (it != points.end())
	{
		index = std::distance(points.begin(), it); // Element index
	}
	return index;
}
//Return adjacent nodes around the given node
void CATPGeoModel::GetAdjacentNodes(vtkUnstructuredGrid* grid, vtkIdType nodeId, std::vector<int> &adjacentNodes)
{
	// Find node index
	//vtkIdType nodeIndex = grid->FindPoint(nodeId);

	if (nodeId >= 0)
	{
		// Get all cells incident to this node
		vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
		grid->GetPointCells(nodeId, cellIds);

		// Iterate each cell
		for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
		{
			vtkIdType cellId = cellIds->GetId(i);

			// Get point ids of the current cell
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
			grid->GetCellPoints(cellId, cellPointIds);

			// Add point ids to adjacent list
			for (vtkIdType j = 0; j < cellPointIds->GetNumberOfIds(); j++)
			{
				vtkIdType pointId = cellPointIds->GetId(j);

				// Skip the current node
				if (pointId != nodeId && std::find(adjacentNodes.begin(), adjacentNodes.end(), pointId) == adjacentNodes.end())
				{
					adjacentNodes.emplace_back(pointId);
				}
			}
		}
	}
}
vtkSmartPointer<vtkUnstructuredGrid> CATPGeoModel::Floater(vtkUnstructuredGrid* grid, vector<intersetionPoint> keyPoints, WEIGHT_TYPE type)
{
	//Extract the outer boundary polyline of the mesh
	vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
	geometryFilter->SetInputData(grid);
	geometryFilter->PassThroughPointIdsOn();
	geometryFilter->Update();

	vtkPolyData* boundaryPolyData = geometryFilter->GetOutput();
	vtkSmartPointer<vtkFeatureEdges> featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
	featureEdges->SetInputData(boundaryPolyData);

	featureEdges->BoundaryEdgesOn();
	featureEdges->ManifoldEdgesOff();
	featureEdges->NonManifoldEdgesOff();
	featureEdges->FeatureEdgesOff();
	featureEdges->Update();

	//Get boundary nodes. Boundary nodes and fracture key points are two special point sets.
	//Boundary node coordinates remain unchanged (same as the reference surface).
	//Fracture key points are mapped to the corresponding locations on the reference surface.
	std::vector<Point> boundPoints;
	for (int i = 0; i < featureEdges->GetOutput()->GetPoints()->GetNumberOfPoints(); i++)
	{
		double point[3];
		featureEdges->GetOutput()->GetPoints()->GetPoint(i, point);
		boundPoints.emplace_back(Point(point, i));
	}

	//Get adjacent node indices for every node (used for parameterization)
	std::vector<std::vector<int>> vecNodesAdjacent;
	for (int i = 0; i < grid->GetNumberOfPoints(); ++i)
	{
		std::vector<int> vecTemp;
		GetAdjacentNodes(grid, i, vecTemp);
		vecNodesAdjacent.emplace_back(vecTemp);
	}

	//Treat boundary points as known and solve a linear system for interior coordinates
	int totalNum = grid->GetNumberOfPoints();
	SpMat A(totalNum, totalNum);
	// fill A
	Rsv b1(totalNum), b2(totalNum);
	for (int i = 0; i < totalNum; ++i)
	{
		b1(i) = 0;
		b2(i) = 0;
	}

	for (int i = 0; i < totalNum; ++i)
	{
		A.coeffRef(i, i) = 1;//Diagonal entry for node i is always 1

		std::vector<double > iWeight = Weigth(grid->GetPoints(), i, vecNodesAdjacent[i], type);

		int index = findPoint(boundPoints, Point(grid->GetPoint(i)));

		if (index >= 0)//Node i is on the boundary; its index in the ordered boundary list is "index"
		{
			//Set RHS to the known coordinates
			b1[i] = grid->GetPoint(i)[0];
			b2[i] = grid->GetPoint(i)[1];
		}
		else
		{
			auto index = findIndexInPoint(i, keyPoints);
			if (index != -1)//Key point on a fracture: coordinates are known but should map to the reference fracture location
			{
				b1[i] = keyPoints[index].p[0];
				b2[i] = keyPoints[index].p[1];
			}
			else//Regular interior node
			{
				for (int j = 0; j < vecNodesAdjacent[i].size(); ++j)
				{
					int index = findPoint(boundPoints, Point(grid->GetPoint(vecNodesAdjacent[i][j])));

					if (index >= 0)//Adjacent node is on the boundary: its UV coordinate is known; move to RHS
					{
						b1[i] += iWeight[j] * grid->GetPoint(vecNodesAdjacent[i][j])[0];
						b2[i] += iWeight[j] * grid->GetPoint(vecNodesAdjacent[i][j])[1];
					}
					else
					{
						A.coeffRef(i, vecNodesAdjacent[i][j]) += -iWeight[j];
					}
				}
			}
		}
	}
	//Eigen::PardisoLU<Eigen::SparseMatrix<double>> solver;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	if (solver.info() != Eigen::Success) {
		// 
		std::cout << "decomposition failed\n";
		return nullptr;
	}
	Rsv x = solver.solve(b1);
	if (solver.info() != Eigen::Success) {
		// solving failed
		std::cout << "solving failed\n";
		return nullptr;
	}
	// solve for another right hand side:
	Rsv y = solver.solve(b2);
	//x and y are the parameterized node coordinates (for all nodes)
	//Create a VTK grid to store the parameterized mesh
	vtkSmartPointer<vtkUnstructuredGrid> ParameterizationGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	ParameterizationGrid->DeepCopy(grid);
	for (int i = 0; i < ParameterizationGrid->GetNumberOfPoints(); ++i)
	{
		ParameterizationGrid->GetPoints()->SetPoint(i, x[i], y[i], 0);
	}
	return ParameterizationGrid;
}
int CATPGeoModel::GetPointNumberOfFracs()
{
	return m_nPointsNumberOfFracs;
}
Point CATPGeoModel::calculateIntersection(Point A, Point B, Point C, Point D)
{
	struct Line3D {
		double x0, y0, z0; // A point on the line
		double a, b, c;    // Direction vector
	};

	Line3D line1, line2;
	line1.x0 = A._x[0];
	line1.y0 = A._x[1];
	line1.z0 = A._x[2];
	line1.a = B._x[0] - A._x[0];
	line1.b = B._x[1] - A._x[1];
	line1.c = B._x[2] - A._x[2];

	line2.x0 = C._x[0];
	line2.y0 = C._x[1];
	line2.z0 = C._x[2];
	line2.a = D._x[0] - C._x[0];
	line2.b = D._x[1] - C._x[1];
	line2.c = D._x[2] - C._x[2];

	double determinant = line1.a * line2.b - line2.a * line1.b;
	if (0==determinant)
	{
		return { NAN, NAN, NAN };
	}

	double t = ((C._x[0] - A._x[0]) * line2.b - (C._x[1] - A._x[1]) * line2.a) / (line1.a * line2.b - line1.b * line2.a);
	double intersectionX = A._x[0] + t * (B._x[0] - A._x[0]);
	double intersectionY = A._x[1] + t * (B._x[1] - A._x[1]);
	double intersectionZ = A._x[2] + t * (B._x[2] - A._x[2]);

	Point intersectionPoint = { intersectionX, intersectionY, intersectionZ };
	return intersectionPoint;
}
bool CATPGeoModel::isInside(Point A, Point B, Point C, Point D, Point p)
{
	if (p._x[0] >= std::min({ A._x[0],B._x[0],C._x[0], D._x[0] }) &&
		p._x[1] <= std::max({ A._x[1],B._x[1],C._x[1], D._x[1] }))
	{
		return true;
	}
	return false;
}

double CATPGeoModel::distanceTwoPoints(Point p1, Point p2)
{
	return sqrt(pow(p1._x[0]-p2._x[0],2)+ pow(p1._x[1] - p2._x[1], 2));
}

vtkMultiBlockDataSet * CATPGeoModel::fracsGridUpdate(vtkMultiBlockDataSet * layerSurfs)
{
	for (int i = 0; i < m_vecFracs.size(); ++i)
	{
		vtkSmartPointer<vtkUnstructuredGrid> fgrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		vtkSmartPointer<vtkPoints> FracPoints = vtkSmartPointer<vtkPoints>::New();
		for (int j = 0; j < layerSurfs->GetNumberOfBlocks(); ++j)
		{
			auto grid = vtkUnstructuredGrid::SafeDownCast(layerSurfs->GetBlock(j));
			for (int k = 0; k < m_vecFracs[i].vecNodes.size(); ++k)
			{
				auto p = grid->GetPoint(m_vecFracs[i].vecNodes[k]);
				//if (!IsPointAlreadyInserted(FracPoints, p))
				FracPoints->InsertNextPoint(p);
			}

		}
		fgrid->SetPoints(FracPoints);

		//Build cells
		for (int j = 0; j < layerSurfs->GetNumberOfBlocks() - 1; ++j)
		{
			for (int k = 0; k < m_vecFracs[i].vecNodes.size() - 1; ++k)
			{
				vtkIdType quad[4] =
				{
					k + j * m_vecFracs[i].vecNodes.size(),
					k + 1 + j * m_vecFracs[i].vecNodes.size(),
					k + 1 + (j + 1) * m_vecFracs[i].vecNodes.size(),
					k + (j + 1) * m_vecFracs[i].vecNodes.size()
				};
				fgrid->InsertNextCell(VTK_QUAD, 4, quad);
			}
		}
		m_fracsGrid->SetBlock(i, fgrid);
	}
	return m_fracsGrid;
}

void CATPGeoModel::Export2VtkFile(std::string fname, vtkUnstructuredGrid* dataset)
{
	vtkNew<vtkUnstructuredGridWriter> writer;
	writer->SetFileName(fname.c_str());
	writer->SetFileVersion(42);
	writer->SetInputData(dataset);
	writer->Write();
}
void CATPGeoModel::Export2VtmFile(std::string fname, vtkMultiBlockDataSet* dataset)
{
	vtkNew<vtkXMLMultiBlockDataWriter> writer;
	writer->SetFileName(fname.c_str());
	writer->SetInputData(dataset);
	writer->Write();
}

void CATPGeoModel::ExportVecFracs(std::string fname)
{
	ofstream f(fname);
	f << m_vecFracs.size() << "\n";
	for (int i=0;i<m_vecFracs.size();++i)
	{
		f << m_vecFracs[i].p1 << "\t" << m_vecFracs[i].p2 << "\t"<< m_vecFracs[i].beta<<"\n";//Start and end points
		f << m_vecFracs[i].vecFracIntersect.size() << "\n";
		for (int j=0;j< m_vecFracs[i].vecFracIntersect.size();++j)
		{
			f << m_vecFracs[i].vecFracIntersect[j].pI << "\t" << m_vecFracs[i].vecFracIntersect[j].index << "\n";
		}
		f << m_vecFracs[i].vecNodes.size() << "\n";
		for (int j = 0; j < m_vecFracs[i].vecNodes.size(); ++j)
		{
			f << m_vecFracs[i].vecNodes[j] << "\n";
		}
	}
}

void CATPGeoModel::ImportVecFracs(std::string fname)
{
	std::ifstream f(fname);
	if (f.is_open()) { // Check whether the file opened successfully
		std::string line;
		int num = 0;
		std::getline(f, line);
		m_vecFracs.resize(std::stoi(line));
		for (int i=0;i< m_vecFracs.size();++i)
		{
			std::getline(f, line);
			std::vector<std::string> tokens; // Store split fields
			std::stringstream ss(line);
			std::string token;
			char delimiter = '\t'; // Tab delimiter
			while (std::getline(ss, token, delimiter)) {
				tokens.push_back(token); // Push field into vector
			}
			m_vecFracs[i].p1 = std::stoi(tokens[0]); m_vecFracs[i].p2 = std::stoi(tokens[1]); m_vecFracs[i].beta = std::stod(tokens[2]);

			std::getline(f, line);
			m_vecFracs[i].vecFracIntersect.resize(std::stoi(line));
			for (int j=0;j< m_vecFracs[i].vecFracIntersect.size();++j)
			{
				std::getline(f, line);
				std::vector<std::string> tokens; // Store split fields
				std::stringstream ss(line);
				std::string token;
				char delimiter = '\t'; // Tab delimiter
				while (std::getline(ss, token, delimiter)) {
					tokens.push_back(token); // Push field into vector
				}
				m_vecFracs[i].vecFracIntersect[j].pI = std::stoi(tokens[0]);
				m_vecFracs[i].vecFracIntersect[j].index = std::stoi(tokens[1]);
			}
			std::getline(f, line);
			m_vecFracs[i].vecNodes.resize(std::stoi(line));
			for (int j=0;j< m_vecFracs[i].vecNodes.size();++j)
			{
				std::getline(f, line);
				m_vecFracs[i].vecNodes[j] = std::stoi(line);
			}
		}		
	}
	
}

