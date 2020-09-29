void Tri_Mesh::saperateFace()
{
	Allmeshs.clear();
	std::vector <VHandle> vhandles;
	std::vector <FHandle> fhandles;
	for (VIter v = vertices_begin(); v != vertices_end(); ++v)
		vhandles.push_back(v.handle());
	for (FIter f = faces_begin(); f != faces_end(); ++f)
		fhandles.push_back(f.handle());

	// twoBoundEdge pair<pair<faceVecID, faceID>, pair<vid1, vid2>>
	typedef std::pair <std::pair<int, int>, std::pair<int, int>> twoBoundStruct;
	std::vector <twoBoundStruct> twoBoundVec;

	//std::cout << "------------face_saperate------------\n";
	for (int i = 0; i < faceVec.size(); i++) 
	{
		// each Big Face
		Tri_Mesh newMesh;
		std::vector <VHandle> newVH;
		std::map <int, int> checkUse;
		//std::cout << "------face " << i << ":------\n";
		for (int j = 0; j < faceVec[i].size(); j++) 
		{
			FaceHandle temp = fhandles[faceVec[i][j]];
			std::vector <VHandle> newFH;
			for (FVIter fv = fv_iter(temp); fv; ++fv) 
			{
				Point p = point(fv.handle());
				//std::cout << "vertex " << i << ":   " << p << "\n";
				if (checkUse.find(fv.handle().idx()) == checkUse.end())
				{
					VertexHandle tempv = newMesh.add_vertex(Point(p));
					newVH.push_back(tempv);
					newFH.push_back(newVH[newVH.size() - 1]);
					checkUse[fv.handle().idx()] = newVH.size() - 1;
				}
				else {
					newFH.push_back(newVH[checkUse[fv.handle().idx()]]);
				}
			}
			newMesh.add_face(newFH);
		}
		int count = 0;
		newMesh.centerPoint = Point(0, 0, 0);
		std::map <int, bool> useCenter;
		for (FIter f = newMesh.faces_begin(); f != newMesh.faces_end(); ++f)
		{
			int bcount = 0;
			Point temp = Point(0, 0, 0);
			/*for (FVIter fv = newMesh.fv_begin(f); fv; ++fv)
			{
				if (newMesh.is_boundary(fv))
					bcount++;
				else temp = newMesh.point(fv);
			}*/

			/*if (bcount == 2)
			{
				newMesh.centerPoint += temp;
				count++;
			}*/
			for (FEIter fe = newMesh.fe_begin(f); fe; ++fe)
			{
				if (newMesh.is_boundary(fe))
				{
					HalfedgeHandle he = newMesh.halfedge_handle(newMesh.edge_handle(fe), 1);
					VHandle ev[2];
					ev[0] = newMesh.to_vertex_handle(he);
					ev[1] = newMesh.from_vertex_handle(he);
					if (useCenter.find(ev[0].idx()) == useCenter.end())
					{
						useCenter[ev[0].idx()] = true;
						newMesh.centerPoint += newMesh.point(ev[0]);
					}
					if (useCenter.find(ev[1].idx()) == useCenter.end())
					{
						useCenter[ev[1].idx()] = true;
						newMesh.centerPoint += newMesh.point(ev[1]);
					}
				}
			}

		}
		//newMesh.centerPoint /= count;
		newMesh.centerPoint /= useCenter.size();


		// new Mesh first processing
		// addMesh: new Mesh first processing
		Tri_Mesh addMesh;
		newVH.clear();
		checkUse.clear();
		newVH.push_back(addMesh.add_vertex(Point(newMesh.centerPoint)));

		//for (FIter f = newMesh.faces_begin(); f != newMesh.faces_end(); ++f) 
		//{
		//	int Boundcount = 0;

		//	std::vector <VHandle> newFH;

		//	for (FVIter fv = newMesh.fv_iter(f); fv; ++fv)
		//	{
		//		if (newMesh.is_boundary(fv))
		//			Boundcount++;
		//	}
		//	if (Boundcount >= 2) {
		//		for (FVIter fv = newMesh.fv_iter(f); fv; ++fv)
		//		{
		//			if (newMesh.is_boundary(fv)) {
		//				if (checkUse.find(fv.handle().idx()) == checkUse.end()) {
		//					newVH.push_back(addMesh.add_vertex(newMesh.point(fv.handle())));
		//					checkUse[fv.handle().idx()] = newVH.size() - 1;
		//					newFH.push_back(newVH[newVH.size() - 1]);
		//				}
		//				else newFH.push_back(newVH[checkUse[fv.handle().idx()]]);
		//			}
		//			else newFH.push_back(newVH[0]);
		//		}
		//		addMesh.add_face(newFH);
		//		
		//	}
		//}

		for (FIter f = newMesh.faces_begin(); f != newMesh.faces_end(); ++f)
		{
			int BoundCount = 0;
			std::vector <VHandle> newFH;
			std::vector <int> boundEdge;
			std::map <int, bool> useVertex;
			for (FEIter fe = newMesh.fe_iter(f); fe; ++fe)
			{
				HalfedgeHandle he = newMesh.halfedge_handle(fe, 0);
				
				useVertex[newMesh.to_vertex_handle(he).idx()];
				useVertex[newMesh.from_vertex_handle(he).idx()];
				if (newMesh.is_boundary(fe))
				{
					BoundCount++;
					boundEdge.push_back(fe.handle().idx());
				}
			}
			if (BoundCount == 1)
			{
				HalfedgeHandle he = newMesh.halfedge_handle(newMesh.edge_handle(boundEdge[0]), 1);
				VHandle ev[2];
				ev[0] = newMesh.to_vertex_handle(he);
				ev[1] = newMesh.from_vertex_handle(he);

				for (FVIter fv = newMesh.fv_iter(f); fv; ++fv) {
					if (fv.handle().idx() == ev[0].idx() || fv.handle().idx() == ev[1].idx()) {
						if (checkUse.find(fv.handle().idx()) == checkUse.end()) {
							newVH.push_back(addMesh.add_vertex(newMesh.point(fv.handle())));
							checkUse[fv.handle().idx()] = newVH.size() - 1;
							newFH.push_back(newVH[newVH.size() - 1]);
						}
						else newFH.push_back(newVH[checkUse[fv.handle().idx()]]);
					}
					else newFH.push_back(newVH[0]);
				}
				addMesh.add_face(newFH);
			}
			else if (BoundCount == 2)
			{
				HalfedgeHandle he0 = newMesh.halfedge_handle(newMesh.edge_handle(boundEdge[0]), 0);
				HalfedgeHandle he1 = newMesh.halfedge_handle(newMesh.edge_handle(boundEdge[1]), 0);
				VHandle ev0[2];
				ev0[0] = newMesh.to_vertex_handle(he0);
				ev0[1] = newMesh.from_vertex_handle(he0);
				VHandle ev1[2];
				ev1[0] = newMesh.to_vertex_handle(he1);
				ev1[1] = newMesh.from_vertex_handle(he1);
				for (FVIter fv = newMesh.fv_iter(f); fv; ) {
					FVIter prefv = fv;
					++fv;
					int preid = prefv.handle().idx();
					int curid = fv.handle().idx();
					bool isEnd = false;
					if (!fv) {
						fv = newMesh.fv_iter(f);
						isEnd = true;
					}

					// test Two BoundEdge
					if (checkUse.find(prefv.handle().idx()) == checkUse.end()) {
						newVH.push_back(addMesh.add_vertex(newMesh.point(prefv.handle())));
						checkUse[prefv.handle().idx()] = newVH.size() - 1;
						newFH.push_back(newVH[newVH.size() - 1]);
					}
					else newFH.push_back(newVH[checkUse[prefv.handle().idx()]]);
					if ((fv.handle().idx() == ev0[0].idx() || fv.handle().idx() == ev0[1].idx())
						&& (prefv.handle().idx() == ev0[0].idx() || prefv.handle().idx() == ev0[1].idx()))
					{

					}
					else if ((fv.handle().idx() == ev1[0].idx() || fv.handle().idx() == ev1[1].idx())
						&& (prefv.handle().idx() == ev1[0].idx() || prefv.handle().idx() == ev1[1].idx()))
					{

					}
					else 
					{
						twoBoundStruct twoBoundData;
						twoBoundData.first.first = i;
						twoBoundData.second.first = checkUse[prefv.handle().idx()];
						twoBoundData.second.second = checkUse[fv.handle().idx()];
						twoBoundVec.push_back(twoBoundData);
					}

					/*if ((fv.handle().idx() == ev0[0].idx() || fv.handle().idx() == ev0[1].idx())
						&& (prefv.handle().idx() == ev0[0].idx() || prefv.handle().idx() == ev0[1].idx())) 
					{
						if (checkUse.find(prefv.handle().idx()) == checkUse.end()) {
							newVH.push_back(addMesh.add_vertex(newMesh.point(prefv.handle())));
							checkUse[prefv.handle().idx()] = newVH.size() - 1;
							newFH.push_back(newVH[newVH.size() - 1]);
						}
						else newFH.push_back(newVH[checkUse[prefv.handle().idx()]]);
						if (checkUse.find(fv.handle().idx()) == checkUse.end()) {
							newVH.push_back(addMesh.add_vertex(newMesh.point(fv.handle())));
							checkUse[fv.handle().idx()] = newVH.size() - 1;
							newFH.push_back(newVH[newVH.size() - 1]);
						}
						else newFH.push_back(newVH[checkUse[fv.handle().idx()]]);
						newFH.push_back(newVH[0]);
						addMesh.add_face(newFH);
						newFH.clear();
					}
					else if ((fv.handle().idx() == ev1[0].idx() || fv.handle().idx() == ev1[1].idx())
						&& (prefv.handle().idx() == ev1[0].idx() || prefv.handle().idx() == ev1[1].idx()))
					{
						if (checkUse.find(prefv.handle().idx()) == checkUse.end()) {
							newVH.push_back(addMesh.add_vertex(newMesh.point(prefv.handle())));
							checkUse[prefv.handle().idx()] = newVH.size() - 1;
							newFH.push_back(newVH[newVH.size() - 1]);
						}
						else newFH.push_back(newVH[checkUse[prefv.handle().idx()]]);
						if (checkUse.find(fv.handle().idx()) == checkUse.end()) {
							newVH.push_back(addMesh.add_vertex(newMesh.point(fv.handle())));
							checkUse[fv.handle().idx()] = newVH.size() - 1;
							newFH.push_back(newVH[newVH.size() - 1]);
						}
						else newFH.push_back(newVH[checkUse[fv.handle().idx()]]);
						newFH.push_back(newVH[0]);
						addMesh.add_face(newFH);
						newFH.clear();
					}*/
					if (isEnd) break;
				}
				FaceHandle fh = addMesh.add_face(newFH);
				twoBoundVec[twoBoundVec.size() - 1].first.second = fh.idx();
			}
			else if (BoundCount == 3)
			{
				for (FVIter fv = newMesh.fv_iter(f); fv; ++fv)
				{
					if (newMesh.is_boundary(fv)) {
						if (checkUse.find(fv.handle().idx()) == checkUse.end()) {
							newVH.push_back(addMesh.add_vertex(newMesh.point(fv.handle())));
							checkUse[fv.handle().idx()] = newVH.size() - 1;
							newFH.push_back(newVH[newVH.size() - 1]);
						}
						else newFH.push_back(newVH[checkUse[fv.handle().idx()]]);
					}
					else newFH.push_back(newVH[0]);
				}
				addMesh.add_face(newFH);
			}
		}

		// new Mesh second processing
		std::map <int, std::pair<int, std::map<int, bool>>> boundaryConnect;
		for (VIter v = addMesh.vertices_begin(); v != addMesh.vertices_end(); ++v)
		{
			if (addMesh.is_boundary(v))
			{
				boundaryConnect[v.handle().idx()] = std::pair<int, std::map<int, bool>>();
				boundaryConnect[v.handle().idx()].first = v.handle().idx();
				for (VVIter vv = addMesh.vv_begin(v); vv; ++vv)
				{
					if (addMesh.is_boundary(vv))
					{
						boundaryConnect[v.handle().idx()].second[vv.handle().idx()] = true;
					}
				}
			}
		}
		//if (boundaryConnect.find(connect[i]) == boundaryConnect.end())
		//{
		//	boundaryConnect[connect[i]] = std::pair<int, std::map<int, bool>>();
		//	boundaryConnect[connect[i]].first = connect[i];
		//}
		//if (i != j)
		//{
		//	boundaryConnect[connect[i]].second[connect[j]] = true;
		//}

		Allmeshs.push_back(addMesh);
	}
	int count = 0;
	for (int i = 0; i < twoBoundVec.size(); i++)
	{
		int faceVecID = twoBoundVec[i].first.first;
		int faceID = twoBoundVec[i].first.second;
		int v1 = twoBoundVec[i].second.first;
		int v2 = twoBoundVec[i].second.second;
		Tri_Mesh* temp = &Allmeshs[faceVecID];
		if (faceID == -1) continue;
		for (FEIter fe = temp->fe_begin(temp->face_handle(faceID)); fe; ++fe)
		{
			int feid = fe.handle().idx();
			HalfedgeHandle he = temp->halfedge_handle(fe, 0);
			VHandle ev[2];
			ev[0] = temp->to_vertex_handle(he);
			ev[1] = temp->from_vertex_handle(he);
			if ((ev[0].idx() == v1 || ev[0].idx() == v2)
				&& (ev[1].idx() == v1 || ev[1].idx() == v2))
			{
				FaceHandle f = temp->opposite_face_handle(he);
				int fid = f.idx();
				if (fid == -1 || fid == faceID)
				{
					std::vector<VHandle> newFH;
					newFH.push_back(temp->vertex_handle(v2));
					newFH.push_back(temp->vertex_handle(v1));
					newFH.push_back(temp->vertex_handle(0));
					temp->add_face(newFH);
				}
			}
		}
	}
	//std::cout << "------------face_saperate------------\n";
}