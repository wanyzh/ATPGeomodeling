#pragma once
#pragma warning(disable: 4996) // Disable compiler warning C4996
#include <iostream>
#include <vector>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSet.h>
using namespace std;
#define _USE_MATH_DEFINES
#include <math.h>
enum WEIGHT_TYPE {
	Tutte,
	Floater_Shape_Preserve,
	Floater_Mean_Value,
	Discrete_Authalic_Parameterization,
	Discrete_Conformal_Map,
	Iterative_Authalic_Parameterization
};

class Point {
public:
	Point(double x, double y, double z) {
		_x[0] = x;
		_x[1] = y;
		_x[2] = z;
	};
	Point(double x, double y, double z, int index) {
		_x[0] = x;
		_x[1] = y;
		_x[2] = z;
		_index = index;
	};
	Point(double x[3])
	{
		for (int i = 0; i < 3; ++i)
		{
			_x[i] = x[i];
		}
	};
	Point(double x[3], int index)
	{
		for (int i = 0; i < 3; ++i)
		{
			_x[i] = x[i];
		}
		_index = index;
	};
	Point() {};
	Point(const Point& obj)
	{
		_index = obj._index;
		for (int i = 0; i < 3; ++i)
		{
			_x[i] = obj._x[i];
		}
	}

	Point GetXRot(double beta);//Compute coordinates after rotating p0 around X by beta
	Point GetZRot(double alpha);//Compute coordinates after rotating p0 around Z by alpha
	Point Move(Point dL);//Compute coordinates after translating p0 by dL

	Point& operator=(const Point& obj)
	{
		if (this != &obj)
		{
			_index = obj._index;
			for (int i = 0; i < 3; ++i)
			{
				_x[i] = obj._x[i];
			}
		}
		return *this;
	}
	bool operator==(const Point& other)//Points are equal if coordinates match; ignore _index
	{
		for (int i = 0; i < 3; ++i)
		{
			if (_x[i] != other._x[i])
			{
				return false;
			}
		}
		return true;
	}
	// Overload subtraction operator
	Point operator-(const Point& other) const {
		double newX[3];
		for (int i = 0; i < 3; ++i)
		{
			newX[i] = _x[i] - other._x[i];
		}
		return Point(newX);
	}
	// Overload multiplication operator for dot product
	double operator*(const Point& other) const {
		return _x[0] * other._x[0] + _x[1] * other._x[1] + _x[2] * other._x[2];
	}
	// Length
	double GetLength()
	{
		return sqrt(pow(_x[0], 2) + pow(_x[1], 2) + pow(_x[2], 2));
	}
	// Set coordinates
	void SetCoor(double x[3])
	{
		for (int i = 0; i < 3; ++i)
		{
			_x[i] = x[i];
		}
	}
	void SetCoor(double x, double y, double z = 0)
	{
		_x[0] = x;
		_x[1] = y;
		_x[2] = z;
	}
	void SetCoor(const Point &p)
	{
		for (int i = 0; i < 3; ++i)
		{
			_x[i] = p._x[i];
		}
	}


public:
	int _index{ -1 };
	double _x[3]{ 0,0,0 };
};

struct intersetionPoint
{
	double p[3];
	int pI;//Index of this node in the mesh for p[3]
	int index;//Fracture index associated with this intersection
	intersetionPoint()
	{

	}
	intersetionPoint(const intersetionPoint &obj)
	{
		for (int i=0;i<3;++i)
		{
			p[i] = obj.p[i];
		}
		pI = obj.pI;
		index = obj.index;
	}
	intersetionPoint(double temP[3])
	{
		for (int i = 0; i < 3; ++i)
		{
			p[i] = temP[i];
		}
	}
	intersetionPoint(int tempI, int tempIndex)
	{		
		pI = tempI;
		index = tempIndex;
	}
	intersetionPoint(double temP[3], int tempI, int tempIndex=-1)
	{
		for (int i=0;i<3;++i)
		{
			p[i] = temP[i];
		}
		pI = tempI;
		index = tempIndex;
	}
	intersetionPoint(Point po, int tempI, int tempIndex = -1)
	{
		for (int i = 0; i < 3; ++i)
		{
			p[i] = po._x[i];
		}
		pI = tempI;
		index = tempIndex;
	}
};

struct bound
{
	vector<int> _p;//Boundary node indices
};

struct frac
{
	int p1;//Endpoint 1
	int p2;//Endpoint 2
	double beta;//Dip angle in radians

	double coP1[3];//Coordinates of the first point on the fracture segment
	double coP2[3];//Coordinates of the second point on the fracture segment

	int index;//When used as a ray: which endpoint of the fracture this ray belongs to
	int frac_index;//When used as a ray: which fracture this ray originates from

	//Intersections with other fractures
	//For intersections: intersetionPoint::p[3] is the intersection coordinate, and index is the other fracture's index
	//Intersections are ordered from p1 to p2
	vector<intersetionPoint> vecFracIntersect;

	//Whether the fracture endpoints lie on the boundary
	//Only endpoints can intersect the boundary; mark the two endpoints
	//0 = not on boundary, 1 = on boundary
	int isOnBound[2]{0,0};

	//All mesh node indices on the fracture (including all intersections), ordered from p1 to p2
	vector<int> vecNodes;

	//All point coordinates on the fracture, aligned with vecNodes
	vector<Point> vecPoints;

	frac()	{}
	frac(const frac &obj)
	{
		p1 = obj.p1; p2 = obj.p2; beta = obj.beta;
		coP1[0] = obj.coP1[0]; coP1[1] = obj.coP1[1]; coP1[2] = obj.coP1[2];
		coP2[0] = obj.coP2[0]; coP2[1] = obj.coP2[1]; coP2[2] = obj.coP2[2];
		index = obj.index;
		frac_index = obj.frac_index;
		vecFracIntersect = obj.vecFracIntersect;
		isOnBound[0] = obj.isOnBound[0]; isOnBound[1] = obj.isOnBound[1];
		vecNodes = obj.vecNodes;
		vecPoints = obj.vecPoints;
	}
};
class CATPGeoModel
{
public:
	CATPGeoModel();
	void SetInputData(vtkMultiBlockDataSet* allSurfs, vector<frac> fracs,vtkUnstructuredGrid* refGrid);
	void SetInputData(vtkMultiBlockDataSet* allSurfs, vtkUnstructuredGrid* refGrid);
	vtkSmartPointer<vtkUnstructuredGrid> Generate();//Generate a 3D prism mesh from input surfaces and reference grid
	vtkSmartPointer<vtkMultiBlockDataSet> GetFracs();//Return multiblock dataset for fracture cells
	vtkSmartPointer<vtkMultiBlockDataSet> GetSurfOutPut();//Return all parameterized original surfaces
	vtkSmartPointer<vtkUnstructuredGrid> RefineLayer(vector<int> num);
	vtkSmartPointer<vtkMultiBlockDataSet> GetSurfRefineOutPut();//Return all refined surfaces
	void Export2VtkFile(std::string fname, vtkUnstructuredGrid* dataset);
	void Export2VtmFile(std::string fnmae, vtkMultiBlockDataSet* MBDataSet);
	void ExportVecFracs(std::string fname);
	void ImportVecFracs(std::string fname);
	vtkMultiBlockDataSet* getSurfRefineOutputMB() {
		return m_surfRefineOutputMB;
	}

private:
	double getAngle(Point center, Point p);
	bool comparePoints(Point center, Point a, Point b);
	std::vector<Point> sortPointsByPolarAngle(std::vector<Point>& points, Point center);
	double sumOfVectorElements(const std::vector<double>& vec);
	bool compareDoubleArrays(const double* array1, const double* array2, int size);
	bool compareTwoFrac(frac f1, frac f2);
	bool IsfracInvectors(frac f, vector<frac> fracs);
	bool IsPointAlreadyInserted(vtkPoints* points, double point[3]);
	vtkUnstructuredGrid* lineGridConsructAndOutput(vector<vector<frac>> vecFracs, string filename);
	vtkMultiBlockDataSet* fracGridConsructAndOutput(vector<vector<frac>> vecFracs, string filename);
	std::vector<double> Weigth(vtkPoints* points, vtkIdType node, std::vector<int> nodeAdjacents, WEIGHT_TYPE type);
	int findIndexInPoint(int index, vector<intersetionPoint> keyPoints);
	int findIndexInPoint(int index, vector<Point> points);
	int findPoint(std::vector<Point> points, Point p);
	void GetAdjacentNodes(vtkUnstructuredGrid* grid, vtkIdType nodeId, std::vector<int> &adjacentNodes);
	vtkSmartPointer<vtkUnstructuredGrid> Floater(vtkUnstructuredGrid* grid, vector<intersetionPoint> keyPoints, WEIGHT_TYPE type = Tutte);
	int GetPointNumberOfFracs();
	Point calculateIntersection(Point A, Point B, Point C, Point D);
	bool isInside(Point A, Point B, Point C, Point D, Point p);
	double distanceTwoPoints(Point p1, Point p2);
	vtkMultiBlockDataSet* fracsGridUpdate(vtkMultiBlockDataSet* layerSurf);
	

	
private:
	vtkSmartPointer<vtkMultiBlockDataSet> m_surfInputMB;//All input original surfaces
	vtkSmartPointer<vtkMultiBlockDataSet> m_surfOutputMB;//All parameterized surfaces
	vector<frac> m_vecFracs;
	vtkSmartPointer<vtkMultiBlockDataSet> m_fracsGrid;//All fracture grids
	vtkSmartPointer<vtkMultiBlockDataSet> m_surfRefineOutputMB;//All refined surfaces

	vtkSmartPointer<vtkUnstructuredGrid> m_referenceGrid;//Reference surface mesh used as template
	int m_nPointsNumberOfFracs;//Number of fracture mesh nodes on one surface

	double m_dMaxLength;
};

