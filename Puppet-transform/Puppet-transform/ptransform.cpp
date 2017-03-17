//
//  ptransform.cpp
//  Puppet-transform
//
//  Created by t_wangju on 09/02/2017.
//  Copyright © 2017 t_wangju. All rights reserved.
//

#include "ptransform.hpp"

using namespace psl2t;
using namespace Eigen;

// argtg x = pi/4 * x +c
// (p,i)x(p,j) = k * sin ai;
// (p,i)*(p,j) = k' * cos ai;
// c = pi / 88.35;

PTrans::PTrans(set<Triangle*> its):ts(its){
	reorder();
	pv.clear();
	p2s.clear();
	constraintPoints.clear();
	
	// points collection and seq;
	set<Point*> points;
	for(auto t:ts) for(auto p:t->v) points.insert(p);
	int seq=0;
	for(auto p:points){
		p2s[p]=seq++;
		pv.push_back(p);
	}
	
	
	// edges collection and seq;(between pv[i] & pv[j] there is an edge)
	ev.clear();
	set<pair<int, int>> edges;
	for(auto p:pv) for(auto t:p->linkT) for(auto v:t->v) if(v!=p)
		edges.insert(make_pair(p2s[p], p2s[v]));
	for(auto e:edges) ev.push_back(e);
	
	// construct Topological relationships between triangles;
	gmv.clear();
	tv.clear();
	for(auto e:ev){
		vector<int> tmp;
		gmv.push_back(gma(e, tmp));
		tv.push_back(tmp);
	}
}

MatrixXf PTrans::gma(pair<int, int> edge,vector<int> &c){
	c.clear();
	Point *vi=pv[edge.first];
	Point *vj=pv[edge.second];
	Point *vl=NULL,*vr=NULL;
	for(auto t:vi->linkT){
		auto iter=vj->linkT.find(t);
		if(iter!=vj->linkT.end()){
			for(auto p:t->v){
				if(p!=vi&&p!=vj){
					if(p->left(vi, vj)>0)
						vl=p;
					else
						vr=p;
				}
			}
		}
	}

	MatrixXf g;
	if(vl==NULL||vr==NULL){
		vl=vl==NULL?vr:vl;
		g.resize(6, 2);
		g<<
		vi->x,	 vi->y,
		vi->y,	-vi->x,
		vj->x,	 vj->y,
		vj->y,	-vj->x,
		vl->x,	 vl->y,
		vl->y,	-vl->x;
		c.push_back(p2s[vi]);
		c.push_back(p2s[vj]);
		c.push_back(p2s[vl]);
	}else{
		g.resize(8, 2);
		g<<
		vi->x,	 vi->y,
		vi->y,	-vi->x,
		vj->x,	 vj->y,
		vj->y,	-vj->x,
		vl->x,	 vl->y,
		vl->y,	-vl->x,
		vr->x,	 vr->y,
		vr->y,	-vr->x;
		c.push_back(p2s[vi]);
		c.push_back(p2s[vj]);
		c.push_back(p2s[vl]);
		c.push_back(p2s[vr]);
	}
	MatrixXf gt=g.transpose();
	MatrixXf t=gt*g;
	MatrixXf gi=t.inverse();
	MatrixXf res=gi*gt;
	return res;
}

void PTrans::HMatrixConstruct(){
	int count=0;
	HMatrix.resize(ev.size()*2+constraintPoints.size()*2, pv.size()*2);
	HMatrix.setZero();
	for(int i=0;i<ev.size();i++){
		Matrix2f e(2,2);
		e<<
		pv[ev[i].second]->x-pv[ev[i].first]->x,
		pv[ev[i].second]->y-pv[ev[i].first]->y,
		pv[ev[i].second]->y-pv[ev[i].first]->y,
		pv[ev[i].first]->x-pv[ev[i].second]->x;
		MatrixXf h;
		if(gmv[i].size()==12){
			h.resize(2, 6);
			h<<
			-1,0,1,0,0,0,
			0,-1,0,1,0,0;
		}else if(gmv[i].size()==16){
			h.resize(2, 8);
			h<<
			-1,0,1,0,0,0,0,0,
			0,-1,0,1,0,0,0,0;
		}
		MatrixXf r=h-e*gmv[i];
		for(int j=0;j<tv[i].size();j++){
			HMatrix(count*2+0,tv[i][j]*2+0)=r(0,j*2+0);
			HMatrix(count*2+0,tv[i][j]*2+1)=r(0,j*2+1);
			HMatrix(count*2+1,tv[i][j]*2+0)=r(1,j*2+0);
			HMatrix(count*2+1,tv[i][j]*2+1)=r(1,j*2+1);
		}
		count++;
	}
	for(auto p:constraintPoints){
		HMatrix(count*2+0,p2s[p]*2+0)=weight;
		HMatrix(count*2+1,p2s[p]*2+1)=weight;
		count++;
	}
	MatrixXf HT=HMatrix.transpose();
	MatrixXf HP=HT*HMatrix;
	MatrixXf HI=HP.inverse();
	FMatrix=HI*HT;
}

void PTrans::flush(){
	MatrixXf B(ev.size()*2+constraintPoints.size()*2,1);
	B.setZero();
	for(int i=(int)ev.size();i<ev.size()+constraintPoints.size();i++){
		B(i*2+0,0)=weight*constraintPoints[i-ev.size()]->x;
		B(i*2+1,0)=weight*constraintPoints[i-ev.size()]->y;
	}
	pv_similarity.clear();
	
	MatrixXf res=FMatrix*B;
	
	int count=0;
	HMatrix.resize(ev.size()*2+constraintPoints.size()*2, pv.size()*2);
	HMatrix.setZero();
	B.setZero();
	for(int i=0;i<ev.size();i++){
		HMatrix(count*2+0,ev[i].first*2+0)=-1;
		HMatrix(count*2+0,ev[i].second*2+0)=1;
		HMatrix(count*2+1,ev[i].first*2+1)=-1;
		HMatrix(count*2+1,ev[i].second*2+1)=1;
		float ex=pv[ev[i].second]->x-pv[ev[i].first]->x;
		float ey=pv[ev[i].second]->y-pv[ev[i].first]->y;
		MatrixXf e(2,1);
		e<<ex,ey;
		
		MatrixXf tf(2,1);
		MatrixXf b(tv[i].size()*2,1);
		for(int j=0;j<tv[i].size();j++){
			b(j*2+0,0)=res(tv[i][j]*2+0,0);
			b(j*2+1,0)=res(tv[i][j]*2+1,0);
		}
		tf=gmv[i]*b;
		MatrixXf tr(2,2);
		tr<<
		tf(0,0),	tf(1,0),
		-tf(1,0),	tf(0,0);
		tr=tr/sqrt(tf(0,0)*tf(0,0)+tf(1,0)*tf(1,0));
		e=tr*e;
		B(count*2+0,0)=e(0,0);
		B(count*2+1,0)=e(1,0);
		count++;
	}
	
	for(auto p:constraintPoints){
		HMatrix(count*2+0,p2s[p]*2+0)=weight;
		HMatrix(count*2+1,p2s[p]*2+1)=weight;
		B(count*2+0,0)=weight*p->x;
		B(count*2+1,0)=weight*p->y;
		count++;
	}
	
//	res=((((HMatrix.transpose())*HMatrix).inverse())*(HMatrix.transpose()))*B;
//	res=(HMatrix.transpose()*HMatrix).piv
	pv_similarity.clear();
	for(int i=0;i<pv.size();i++){
		pv_similarity.push_back(Point(res(i*2+0,0),res(i*2+1,0)));
	}
}

void PTrans::reorder(){
	for(auto t:ts){
		if(t->v[2]->left(t->v[0], t->v[1])==-1){
			swap(t->v[0], t->v[1]);
			cout<<endl;
		}
	}
}



