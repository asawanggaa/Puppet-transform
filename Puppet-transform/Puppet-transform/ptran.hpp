//
//  ptran.hpp
//  Puppet-transform
//
//  Created by t_wangju on 21/02/2017.
//  Copyright © 2017 t_wangju. All rights reserved.
//

#ifndef ptran_hpp
#define ptran_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <set>
#include <eigen3/Eigen/Dense>

#include "shape.hpp"

using namespace std;
using namespace Eigen;
using namespace psl2t;

class ptran{
public:
	ptran(set<Triangle*> ts);
	void InsertConstraintPoints(Point *p);
	
	// these group of functions are for inside compute;
//private:
	void SetConstriant();
	void FlushConstriant();
	
	// these group of data are for input and initial;
//private:
	vector<Point*> Points;
	vector<Point*> ConstraintPoints;
	vector<Triangle*> Triangles;
	map<Point*,int> PointsOrder;
	map<Triangle*,int> TrianglesOrder;
	// these group of data are for memorized intermidiate data;
//private:
	VectorXf IntermidiateV;
	VectorXf FittedV;
	
	// these group of data are for flush graphic;
//private:
	// Es=v'Gv;
	MatrixXf PointsTopologyM;// topology relationship between points;
	MatrixXf SimilarityM;// -G'.inverse*B for flush operation;
	vector<Point*> SimilarityPoints;
	
	map<Triangle*,vector<Point>> shadow;
	// Ef=v'Hv+fv+c;
	MatrixXf EdgesTopologyM;// topology relationship between edges;
	MatrixXf FinallyM;
	MatrixXf CongruentM;//
	
	
	// partial derivatives of energy; Fw+C=0, memorized F.inverse;
	vector<MatrixXf> RotateMInverseV;
};

#endif /* ptran_hpp */
