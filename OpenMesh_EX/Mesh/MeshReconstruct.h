#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "FaceData.h"
#include "Clustering.h"



inline void  filtminf(std::vector<faceData> fd, std::vector<Cluster>& output )
{
	double area = 0;

	for (int i = 0; i < output.size(); i++)
	{
		for (int j = 0; j < output[i].faceid.size(); j++)
		{
			//std::cout << " face size " << output[i].faceid.size() << std::endl;
			area = fd[output[i].faceid[j]].getArea();
			if (area < (double)0.6)
			{
				output[i].faceid.erase(output[i].faceid.begin() + j);
			//	std::cout << " face size " << output[i].faceid.size() << std::endl;
				j = 0;
				if (output[i].faceid.size() == 1)
				{
					output[i].faceid.pop_back();
				}
			}
		}
		if (output[i].faceid.size() == 0)
		{
			output.erase(output.begin() + i);
			i = 0;
		}
	}
	if (output.size() == 1)
	{
		output.pop_back();
	}
}

inline void Refine( std::vector<faceData> fd, std::vector<Cluster>& cl )
{
	
	bool merged = false;
	float theta = 10.0f;
	theta  = static_cast<float>(M_PI * theta / 180.0f);
	do 
	{
		merged = false;
		for(int i= 0 ; i<cl.size();i++ )
		{
			for (int j = i + 1; j < cl.size(); ++j)
			{
				if (OpenMesh::dot(fd[cl[i].faceid[0]].getFaceNormal() ,fd[cl[j].faceid[0]].getFaceNormal()) <0.9)
				{
					if (abs(OpenMesh::dot(fd[cl[i].faceid[0]].getFaceNormal(), fd[cl[j].faceid[0]].getFaceNormal())) < std::cos(theta))
					{
						int idx;
					
						Cluster tmp;

							 tmp = cl[idx];
							cl.erase(cl.begin() + idx);

						for (int k = 0; k < tmp.faceid.size(); k++)
						{
							cl[i].faceid.push_back(tmp.faceid[k]);
						}
						merged = true;
						break;
					}
				}
			}
		}
		if (merged)
			break;
	} while (merged);

	
}

inline void hypothesisplane(std::vector<faceData> fd, std::vector<Cluster>& cl)
{

	OpenMesh::Vec3d norm1 = fd[cl[0].faceid[0]].getFaceNormal();
	OpenMesh::Vec3d cent1 = fd[cl[0].faceid[0]].getfcenter();
	double xmax= -10000, xmin = 10000, ymax = -10000, ymin = 10000, zmax = -1000, zmin= 1000;
	/*
	for (int i = 0; i < 1; i++) 
	{
		for (int j = 0; j < cl[i].faceid.size(); j++) 
		{
			for(int k = 0 ; k <3 ; k++)
			{
				if (fd[cl[i].faceid[j]].getVertices()[k][0] >xmax) 
				{
					xmax = fd[cl[i].faceid[j]].getVertices()[k][0];
				}
				if (fd[cl[i].faceid[j]].getVertices()[k][0] < xmin) 
				{
					xmin = fd[cl[i].faceid[j]].getVertices()[k][0];
				}
				//----------------
				if (fd[cl[i].faceid[j]].getVertices()[k][1] > ymax)
				{
					ymax = fd[cl[i].faceid[j]].getVertices()[k][1];
				}
				if (fd[cl[i].faceid[j]].getVertices()[k][1] < ymin)
				{
					ymin = fd[cl[i].faceid[j]].getVertices()[1][0];
				}
				//---------------------
				if (fd[cl[i].faceid[j]].getVertices()[k][2] > zmax)
				{
					zmax = fd[cl[i].faceid[j]].getVertices()[k][2];
				}
				if (fd[cl[i].faceid[j]].getVertices()[k][2] < zmin)
				{
					zmin = fd[cl[i].faceid[j]].getVertices()[k][2];
				}
			}
			
		}
	}
	*/

	OpenMesh::Vec3d lu(xmax, ymax, zmax);
	OpenMesh::Vec3d rd(xmin, ymax, zmin);
	OpenMesh::Vec3d ld(xmax, ymin, zmin);
	OpenMesh::Vec3d ru(xmin, ymin, zmax);

	std::cout << lu << std::endl;
	std::cout << ld << std::endl;
	std::cout << ru << std::endl;
	std::cout << rd << std::endl;

}

