//
//  ptransform.hpp
//  Puppet-transform
//
//  Created by t_wangju on 09/02/2017.
//  Copyright © 2017 t_wangju. All rights reserved.
//

#ifndef ptransform_hpp
#define ptransform_hpp

#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>

#include <eigen3/Eigen/Dense>
#include "shape.hpp"

using namespace psl2t;
using namespace Eigen;



class PTrans{
public:
	PTrans(set<Triangle*> its);

	PTrans(Point* control, bool exist);
	
	void HMatrixConstruct();
	void TopologicalMesh();
	MatrixXf gma(pair<int, int> edge,vector<int> &c);
	void flush();
	
	
	MatrixXf HMatrix;
	MatrixXf FMatrix;
	MatrixXf SMatrix;
	MatrixXf GMatrix;
	
	set<Triangle*> ts;
	void reorder();
	vector<Point*> pv;
	map<Point*,int> p2s;
	vector<Point*> constraintPoints;
	
	vector<pair<int, int>> ev;
	
	vector<vector<int>> tv;
	vector<MatrixXf> gmv;
	
	vector<Point> pv_similarity;
	
	int weight=2000;
};

#endif /* ptransform_hpp */
