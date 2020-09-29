#pragma once
#ifndef _FACE_DATA_
#define _FACE_DATA_
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <utility>
#define pi 3.1415926
static float angleError = 2;
static double normalError = 0.85;
static double faceError = 0.01;
static double faceScale = 1;
static double faceIgnore = 0;

struct faceColor {
	faceColor() {}
	faceColor(float angle, OpenMesh::Vec3d normal) { this->angle = angle; faceNormal = normal; }
	OpenMesh::Vec3d Color;
	int ColorID;
	float angle;
	OpenMesh::Vec3d faceNormal;
	OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DefaultTraits>::FaceIter f_it;
};

bool closeAngle(float a1, float a2);
bool closeNormal(OpenMesh::Vec3d n1, OpenMesh::Vec3d n2);
bool closeFace(std::vector <OpenMesh::Vec3d> vs1, std::vector <OpenMesh::Vec3d> vs2);
bool faceConnect(std::map <int, bool> vid1, std::map<int, bool> vid2);
void setAngleError(float error);
void setNormalError(double error);
void setFaceError(double error);
void setScale(double scale);
double getScale();
void setFaceIgnore(double error);
double getFaceIgnore();

class faceData
{
public:

	faceData();
	faceData(int id);
	~faceData();
	void addvertex(OpenMesh::Vec3d vertex, int vid);
	void caclcenter();
	void setFaceNormal(OpenMesh::Vec3d normal, OpenMesh::Vec3d symmetric);
	void cal_area();
	bool checkData();
	int getID();
	float getAngle();
	double getArea();
	bool isAllow();
	void setAllowFace(bool allow);
	std::vector <OpenMesh::Vec3d> getVertices();
	OpenMesh::Vec3d getFaceNormal();
	OpenMesh::Vec3d getfcenter();
	friend std::ostream& operator<<(std::ostream& os, const faceData& face);
	OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DefaultTraits>::FaceIter f_it;
	std::map <int, bool> vID;
private:
	std::vector <OpenMesh::Vec3d> vertices;
	OpenMesh::Vec3d center;
	OpenMesh::Vec3d faceNormal;
	int id;
	float angle;
	double area;
	bool dataStatus;
	bool allowFace;
};

#endif // !_FACE_COLOR_

