#include "FaceData.h"

faceData::faceData()
{
	id = -1;
	angle = -1;
	area = -1;
	dataStatus = false;
}

faceData::faceData(int id)
{
	this->id = id;
}

faceData::~faceData()
{
	vertices.clear();
	std::vector <OpenMesh::Vec3d>().swap(vertices);
}

void faceData::addvertex(OpenMesh::Vec3d vertex, int vid)
{
	vertices.push_back(vertex);
	vID[vid] = true;
}
void faceData::caclcenter()
{
	if (vertices.size() == 3)
	{
		center[0] = (vertices[0][0] + vertices[1][0] + vertices[2][0]) / 3;
		center[1] = (vertices[0][1] + vertices[1][1] + vertices[2][1]) / 3;
		center[2] = (vertices[0][2] + vertices[1][2] + vertices[2][2]) / 3;
	}
	else 
	{
		center[0] = (vertices[0][0] + vertices[1][0] + vertices[2][0] + vertices[3][0] ) /4;
		center[1] = (vertices[0][1] + vertices[1][1] + vertices[2][1] + vertices[3][1]) / 4;
		center[2] = (vertices[0][2] + vertices[1][2] + vertices[2][2] + vertices[3][2]) / 4;
	}
}
void faceData::setFaceNormal(OpenMesh::Vec3d normal, OpenMesh::Vec3d symmetric)
{
	faceNormal = normal * -1;
	angle = acos(OpenMesh::dot(faceNormal.normalize_cond(), symmetric) / faceNormal.normalize_cond().length() / symmetric.length()) * (180.0f / pi);
	if ((float)angle > 90) angle = 180 - angle;
}

void faceData::cal_area()
{
	double a = (vertices[0] - vertices[1]).length();
	double b = (vertices[1] - vertices[2]).length();
	double c = (vertices[0] - vertices[2]).length();
	double s = (a + b + c) / 2.0;
	area = sqrt(s * (s - a) * (s - b) * (s - c));
	if (area < faceIgnore) allowFace = false;
	else allowFace = true;
}

bool faceData::checkData()
{
	if (id == -1 || angle == -1 || area == -1) {
		dataStatus = false;
		return false;
	}
	dataStatus = true;
	return true;
}

int faceData::getID()
{
	return id;
}

float faceData::getAngle()
{
	return angle;
}

double faceData::getArea()
{
	return area;
}

bool faceData::isAllow()
{
	return allowFace;
}

void faceData::setAllowFace(bool allow)
{
	allowFace = allow;
}

std::vector<OpenMesh::Vec3d> faceData::getVertices()
{
	return vertices;
}

OpenMesh::Vec3d faceData::getFaceNormal()
{
	return faceNormal;
}

OpenMesh::Vec3d faceData::getfcenter()
{
	return center;
}

bool closeAngle(float a1, float a2)
{
	if (abs(a1 - a2) <= angleError) return true;
	return false;
}

bool closeNormal(OpenMesh::Vec3d n1, OpenMesh::Vec3d n2)
{
	if (OpenMesh::dot(n1, n2) >= normalError) return true;
	return false;
}

bool closeFace(std::vector<OpenMesh::Vec3d> vs1, std::vector<OpenMesh::Vec3d> vs2)
{
	for (int i = 0; i < vs1.size(); i++)
		for (int j = 0; j < vs2.size(); j++)
			if ((vs1[i] - vs2[j]).length() <= faceError)
			{
				double len = (vs1[i] - vs2[j]).length();
				return true;
			}
	return false;
}

bool faceConnect(std::map<int, bool> vid1, std::map<int, bool> vid2)
{
	std::map<int, bool>::iterator it1;
	for (it1 = vid1.begin(); it1 != vid1.end(); it1++)
	{
		if (vid2.find(it1->first) != vid2.end())
			return true;
	}
	return false;
}

void setAngleError(float error)
{
	angleError = error;
}

void setNormalError(double error)
{
	normalError = error;
}

void setFaceError(double error)
{
	faceError = error;
}

void setScale(double scale)
{
	faceScale = scale;
}
double getScale()
{
	return faceScale;
}

void setFaceIgnore(double error)
{
	faceIgnore = error;
}

double getFaceIgnore()
{
	return faceError;
}

std::ostream& operator<<(std::ostream& os, const faceData& face)
{
	if (!face.dataStatus)
		os << "Incomplete Data\n";
	else
		os << "Face id: " << face.id << "   face Normal: " << face.faceNormal << "   angle: " << face.angle << "   area: " << face.area;
	os << std::endl;
	return os;
}
