#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "FaceData.h"

typedef struct  Cluster
{
	OpenMesh::Vec3d norm;
	int count = 0;
	std::vector<int> faceid;
	double facearea = 0;
}Cluster;


inline double calecldistance(OpenMesh::Vec3d cur, OpenMesh::Vec3d gl)
{
	// sqrt(pow((gl[0] - cur[0]), 2) + pow((gl[1] - cur[1]), 2) + pow((gl[2] - cur[2]), 2));
	//std::cout << "dist  " << sqrt(pow((gl[0] - cur[0]), 2) + pow((gl[1] - cur[1]), 2) + pow((gl[2] - cur[2]), 2)) << std::endl;
	return pow((gl[0] - cur[0]), 2) + pow((gl[1] - cur[1]), 2) + pow((gl[2] - cur[2]), 2);
}

inline double calplanedist( OpenMesh::Vec3d normal , OpenMesh::Vec3d center,OpenMesh::Vec3d qpoint)
{
	double  planeequation;
	double normalize;
	double dist;

		planeequation = normal[0]*(qpoint[0]-center[0]) + normal[1] * (qpoint[1] - center[1]) + normal[2] * (qpoint[2] - center[2]);
		normalize = sqrtf(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2));
		dist =abs( planeequation / normalize);
		return dist;
}


inline void  kmeans_init(int size, int K, int* pick)
{
	int i, j, k, rnd;
	// 隨機找K 個不同資料
	std::cout << " K :" << K << std::endl;
	for (k = 0; k < K; ++k) {
		//	std::cout<<"idx  " <<k<<"  val  " << pick[k] << std::endl;
		rnd = rand() % size; // 隨機取一筆
		//std::cout << "rnd  " << rnd<<std::endl;
		for (i = 0; i < k && pick[i] != rnd; ++i);
		if (i == k) pick[k] = rnd; // 沒重覆
		else --k; // 有重覆, 再找一次
	}

}


inline double update_table(int* ch_pt, int* cent_c, OpenMesh::Vec3d* dis_k, int K, int len, Cluster t, std::vector<OpenMesh::Vec3d> center, std::vector<faceData>fd, int* table)
{
	int i, j, k, min_k;
	double dis, min_dis, t_sse = 0.0;

	*ch_pt = 0;                          // 變動點數設0
	memset(cent_c, 0, sizeof(cent_c)); // 各叢聚資料數清0
	memset(dis_k, 0, sizeof(dis_k));   // 各叢聚距離和清0

	for (i = 0; i < len; ++i) {
		// 尋找所屬重心
		min_dis = calecldistance(fd[t.faceid[i]].getfcenter(), center[0]);
		min_k = 0;
		for (k = 1; k < K; ++k) {
			dis = calecldistance(fd[t.faceid[i]].getfcenter(), center[k]);
			if (dis < min_dis)
				min_dis = dis, min_k = k;
		}
		*ch_pt += (table[i] != min_k); // 更新變動點數
		table[i] = min_k;          // 更新所屬重心
		++cent_c[min_k];           // 累計重心資料數        
		t_sse += min_dis;          // 累計總重心距離
		for (j = 0; j < 3; ++j)     // 更新各叢聚總距離
			dis_k[min_k][j] += fd[t.faceid[i]].getfcenter()[j];

	}
	return t_sse;
}

inline double caclclusterdis(std::vector<OpenMesh::Vec3d> center)
{
	double dist = 0;

	for (int i = 0; i < center.size(); i++)
	{
		for (int j = i + 1; j < center.size(); j++)
		{
			dist += calecldistance(center[i], center[j]);
		}
	}
	return dist;
}
inline void update_cent(int K, std::vector<OpenMesh::Vec3d>center, int* cent_c, OpenMesh::Vec3d* dis_k)
{
	int k, j;
	for (k = 0; k < K; ++k)
		for (j = 0; j < 3; ++j)
			center[k][j] = dis_k[k][j] / cent_c[k];
}

inline void Kmeans(std::vector<faceData>fd, Cluster* t, int* index, Cluster* output, int* outputidx, std::vector<OpenMesh::Vec3d> ctr)
{
	std::vector<OpenMesh::Vec3d> center;
	int len;
	int Max_iter = 20;
	int Min_PT = 0;
	int* pick;
	int K;
	int* cent_c;
	OpenMesh::Vec3d* dis_k;
	int* table;
	int ch_pt;
	int iter = 0;
	double sse1;
	double sse2;
	double minclustdist = 1000;
	int bk;
	// inti kmeans initial center
	for (int n = 0; n < *index; n++)
	{

		len = t[n].faceid.size();

		if (len > 2) { K = 2; }
		else
		{
			K = 1;
		}

		sse1 = 0;
		sse2 = 0;
		cent_c = new int[K];
		table = new int[len];
		pick = new int[K];
		dis_k = new OpenMesh::Vec3d[K];
		iter = 0;
		center.clear();
		std::cout << "n    " << n << std::endl;
		kmeans_init(len, K, pick);

		for (int k = 0; k < K; k++)
		{

			center.push_back(fd[t[n].faceid[pick[k]]].getfcenter());
		}

		//----------------------------------------
		sse2 = update_table(&ch_pt, cent_c, dis_k, K, len, t[n], center, fd, table);
		do {

			std::cout << "ch_pt  " << ch_pt << std::endl;
			sse1 = sse2, ++iter;
			update_cent(K, center, cent_c, dis_k);
			sse2 = update_table(&ch_pt, cent_c, dis_k, K, len, t[n], center, fd, table);
		} while (iter<Max_iter && sse1 != sse2 && ch_pt>Min_PT);
		std::cout << "testing kmeans" << std::endl;
		std::cout << "sse    =  " << sse2 << std::endl;
		//	std::cout << "ch_pt =  " << ch_pt << std::endl;
			//std::cout << "iter     =  " << iter << std::endl;


		bool isempty = true;
		std::cout << n << "   <-n clusters   =  " << bk << std::endl;
		for (int i = 0; i < K; i++)
		{
			std::cout << "Im inside" << std::endl;
			//t[*index].faceid.push_back(pick[i]);
			for (int j = 0; j < len; j++)
			{
				if (pick[table[j]] == pick[i])
				{
					isempty = false;
					output[*outputidx].faceid.push_back(t[n].faceid[j]);
				}
			}
			if (!isempty)
			{
				*outputidx += 1;
				isempty = true;
			}
		}
	}

	/*
	for (int i = 0; i < 3; i++)
	{
		std::cout << "CLuster  " << pick[i] << "  ------------------------------" << std::endl;
		for (int j = 0; j < len; j++)
		{
			if (pick[table[j]] == pick[i])
			{
				std::cout << j << "    " << pick[table[j]] << std::endl;
			}
		}
	}*/
}



//-----------------------

inline int find_faceid(int targetid , Cluster *cl , int index) 
{
		for (int i = 0; i < index; i++)
		{
			for (int j = 0; j < cl[i].faceid.size(); j++)
			{
				if (cl[i].faceid[j] == targetid) {
					return i;
				}
			}
		}

	return -1;
}

inline int find_normal(Cluster* norcarray, OpenMesh::Vec3d curvec, int size)
{
	for (int i = 0; i < size; i++) {
		if (OpenMesh::dot(curvec, norcarray[i].norm) >= 0.90)
		{
			return i;
		}
	}
	return -1;
}

inline int pcompare(const void* a, const void* b)//這函式是 qsort 所需的比較函式
{
	Cluster va, vb;

	va = *(Cluster*)a;
	vb = *(Cluster*)b;

	if (va.count < vb.count) { return 1; }               //傳回 -1 代表 a < b
	else if (va.count == vb.count) { return 0; }      //傳回   0 代表 a = b
	else return -1;                          //傳回  1 代表 a>b
}

inline int acompare(const void* a, const void* b)//這函式是 qsort 所需的比較函式
{
	Cluster va, vb;

	va = *(Cluster*)a;
	vb = *(Cluster*)b;

	if (va.facearea < vb.facearea) { return 1; }               //傳回 -1 代表 a < b
	else if (va.facearea == vb.facearea) { return 0; }      //傳回   0 代表 a = b
	else return -1;                          //傳回  1 代表 a>b
}