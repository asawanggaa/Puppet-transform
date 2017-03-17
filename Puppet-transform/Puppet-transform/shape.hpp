//
//  shape.hpp
//  Puppet-transform
//
//  Created by t_wangju on 03/02/2017.
//  Copyright © 2017 t_wangju. All rights reserved.
//

#ifndef shape_hpp
#define shape_hpp

#include <stdio.h>
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <queue>
#include <stack>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

namespace psl2t {
	class Point;
	class Triangle;
	class Mesh;
	
	class Point{
	public:
		float x=0.0f;
		float y=0.0f;
		std::set<Triangle*> linkT;
		
		Point(float ix, float iy);
		~Point();
		float dis(const Point &A)const;
		int left(Point *A, Point *B)const;
	};
	
	class Triangle{
	public:
		Point *v[3];
		float coeffi[3][2];
		float R;
		Point O;
		vector<Point> OrdinaryPoints;
		Triangle(Point *A,Point *B,Point *C);
		~Triangle();
		
		bool in(const Point *p)const;
		bool incircle(const Point *p)const;
		Triangle* next(const Point* p)const;
	};
	
	class Mesh{
	public:
		Mesh(float xmin,float ymin,float xmax,float ymax);
		set<Triangle*> Bowyer_Watson(vector<Point*> &input);
//	private:
		set<Triangle*> InsertAPoint(set<Triangle*> &ts, Point *p);
		Triangle* FindFirstBadTriangle(set<Triangle*> &ts,Point* p);
		set<Triangle*> ExploreBadTriangles(set<Triangle*> &ts,Triangle* fbt,Point* p);
		set<set<Point*>> FindUnSharedEdges(set<Triangle*> &badts);
		set<Triangle*> ConstructNewTriangles(set<set<Point*>> &use,Point *p);
		
		set<Triangle*> FindCrossTriangles(set<Triangle*> &ts,Point* s,Point* t);
		
		set<Triangle*> EraseEdgeTriangles();
//	private:
		Point *LB,*RB,*LT,*RT;
		set<Triangle*> ts;
	};
}

#endif /* shape_hpp */
