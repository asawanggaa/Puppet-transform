//
//  ptran.cpp
//  Puppet-transform
//
//  Created by t_wangju on 21/02/2017.
//  Copyright © 2017 t_wangju. All rights reserved.
//

#include "ptran.hpp"
#include <GLUT/GLUT.h>
#include <OpenGL/OpenGL.h>


ptran::ptran(set<Triangle*> ts){
	Points.clear();
	ConstraintPoints.clear();
	Triangles.clear();
	PointsOrder.clear();
	TrianglesOrder.clear();
	
	
	int seq=0;
	set<Point*> ps;
	for(auto t:ts){
		Triangles.push_back(t);
		TrianglesOrder[t]=seq++;
		for(auto v:t->v){
			ps.insert(v);
		}
	}
	
	seq=0;
	for(auto p:ps){
		Points.push_back(p);
		PointsOrder[p]=seq++;
	}
	
	PointsTopologyM.resize(Points.size()*2, Points.size()*2);
	PointsTopologyM.setZero();
	EdgesTopologyM.resize(Points.size()*2, Points.size()*2);
	EdgesTopologyM.setZero();
	for(auto t:Triangles){
		for(int p=0;p<3;p++){
			int i=PointsOrder[t->v[(p+0)%3]];//v0
			int j=PointsOrder[t->v[(p+1)%3]];//v1
			int k=PointsOrder[t->v[(p+2)%3]];//v2
			float x=t->coeffi[(p+2)%3][0];
			float y=t->coeffi[(p+2)%3][1];
			{
				PointsTopologyM(i*2+0,i*2+0)+=1-2*x+x*x+y*y;
				PointsTopologyM(i*2+0,i*2+1)+=0;
				PointsTopologyM(i*2+0,j*2+0)+=2*x-2*x*x-2*y*y;
				PointsTopologyM(i*2+0,j*2+1)+=-2*y;
				PointsTopologyM(i*2+0,k*2+0)+=-2+2*x;
				PointsTopologyM(i*2+0,k*2+1)+=2*y;
				PointsTopologyM(i*2+1,i*2+1)+=1-2*x+x*x+y*y;
				PointsTopologyM(i*2+1,j*2+0)+=2*y;
				PointsTopologyM(i*2+1,j*2+1)+=2*x-2*x*x-2*y*y;
				PointsTopologyM(i*2+1,k*2+0)+=-2*y;
				PointsTopologyM(i*2+1,k*2+1)+=-2+2*x;
				PointsTopologyM(j*2+0,j*2+0)+=x*x+y*y;
				PointsTopologyM(j*2+0,j*2+1)+=0;
				PointsTopologyM(j*2+0,k*2+0)+=-2*x;
				PointsTopologyM(j*2+0,k*2+1)+=-2*y;
				PointsTopologyM(j*2+1,j*2+1)+=x*x+y*y;
				PointsTopologyM(j*2+1,k*2+0)+=2*y;
				PointsTopologyM(j*2+1,k*2+1)+=-2*x;
				PointsTopologyM(k*2+0,k*2+0)+=1;
				PointsTopologyM(k*2+0,k*2+1)+=0;
				PointsTopologyM(k*2+1,k*2+1)+=1;
			}
			{
				EdgesTopologyM(i*2+0,i*2+0)+=1;
				EdgesTopologyM(j*2+0,j*2+0)+=1;
				EdgesTopologyM(i*2+0,j*2+0)+=-2;
				EdgesTopologyM(i*2+1,i*2+1)+=1;
				EdgesTopologyM(j*2+1,j*2+1)+=1;
				EdgesTopologyM(i*2+1,j*2+1)+=-2;
			}
		}
		Matrix4f RM;
		float x=t->coeffi[2][0];
		float y=t->coeffi[2][1];
		RM<<
		2-2*x+x*x+y*y,	0,	x-x*x-y*y,	-y,
		0,	2-2*x+x*x+y*y,	y,	x-x*x-y*y,
		x-x*x-y*y,	y,	1+x*x+y*y,	0,
		-y,	x-x*x-y*y,	0,	1+x*x+y*y;
		Matrix4f RMI=RM.inverse();
		RotateMInverseV.push_back(RMI);
	}
}

void ptran::InsertConstraintPoints(Point *p){
	ConstraintPoints.push_back(p);
}

void ptran::SetConstriant(){
	
	int u=static_cast<int>(Points.size()-ConstraintPoints.size());
	int q=static_cast<int>(ConstraintPoints.size());
	
	vector<int> cons;
	for(auto p:ConstraintPoints) cons.push_back(PointsOrder[p]);
	sort(cons.begin(), cons.end());
	
	auto li=cons.begin();
	auto hi=cons.end()-1;
	int low=0;
	int high=u+q-1;
	
	while(low<high){
		while(low!=*li&&low<high) low++;li++;
		while(high==*hi&&low<high) high--,hi--;
		if(low>=high) break;
		PointsTopologyM.row(low*2+0).swap(PointsTopologyM.row(high*2+0));
		PointsTopologyM.row(low*2+1).swap(PointsTopologyM.row(high*2+1));
		PointsTopologyM.col(low*2+0).swap(PointsTopologyM.col(high*2+0));
		PointsTopologyM.col(low*2+1).swap(PointsTopologyM.col(high*2+1));
		EdgesTopologyM.row(low*2+0).swap(EdgesTopologyM.row(high*2+0));
		EdgesTopologyM.row(low*2+1).swap(EdgesTopologyM.row(high*2+1));
		EdgesTopologyM.col(low*2+0).swap(EdgesTopologyM.col(high*2+0));
		EdgesTopologyM.col(low*2+1).swap(EdgesTopologyM.col(high*2+1));
		swap(PointsOrder[Points[low]], PointsOrder[Points[high]]);
		swap(Points[low], Points[high]);
		high--,hi--;
	}
	MatrixXf GM00(u*2,u*2);
	MatrixXf GM01(u*2,q*2);
	MatrixXf GM10(q*2,u*2);
	GM00=PointsTopologyM.block(0, 0, u*2, u*2);
	GM01=PointsTopologyM.block(0, u*2, u*2, q*2);
	GM10=PointsTopologyM.block(u*2, 0, q*2, u*2);
	MatrixXf GM00T=GM00.transpose();
	MatrixXf G=GM00+GM00T;
	MatrixXf B=GM01+GM10.transpose();
	MatrixXf GI=G.inverse();
	SimilarityM=-(GI*B);
	
	FinallyM.resize(u*2,u*2);
	MatrixXf HM00(u*2,u*2);
	HM00=EdgesTopologyM.block(0, 0, u*2, u*2);
	MatrixXf HM00T=HM00.transpose();
	MatrixXf H=HM00+HM00T;
	FinallyM=H.inverse();
}

void ptran::FlushConstriant(){
	int un=static_cast<int>(Points.size()-ConstraintPoints.size());
	int qn=static_cast<int>(ConstraintPoints.size());
	VectorXf q(qn*2);
	for(int i=0;i<qn;i++){
		q(i*2+0)=Points[un+i]->x;
		q(i*2+1)=Points[un+i]->y;
	}
	VectorXf u=SimilarityM*q;
//	for(int i=0;i<un;i++){
//		Points[i]->x=u(i*2+0);
//		Points[i]->y=u(i*2+1);
//	}
	IntermidiateV.resize(Points.size()*2);
	for(int i=0;i<u.size();i++) IntermidiateV(i)=u(i);
	for(int i=0;i<q.size();i++) IntermidiateV(i+u.size())=q(i);

	shadow.clear();
	for(int i=0;i<Triangles.size();i++){
		Triangle* t=Triangles[i];
		int v0=PointsOrder[t->v[0]];
		int v1=PointsOrder[t->v[1]];
		int v2=PointsOrder[t->v[2]];
		float x=t->coeffi[2][0];
		float y=t->coeffi[2][1];
		Vector4f C;
		C<<
		-IntermidiateV(v0*2+0)-(1-x)*IntermidiateV(v2*2+0)+y*IntermidiateV(v2*2+1),
		-IntermidiateV(v0*2+1)-y*IntermidiateV(v2*2+0)-(1-x)*IntermidiateV(v2*2+1),
		-IntermidiateV(v1*2+0)-x*IntermidiateV(v2*2+0)-y*IntermidiateV(v2*2+1),
		-IntermidiateV(v1*2+1)+y*IntermidiateV(v2*2+0)-x*IntermidiateV(v2*2+1);
		Vector4f w=-RotateMInverseV[i]*C;
		Point p0(w(0),w(1));
		Point p1(w(2),w(3));
		vector<Point> pv;
		pv.push_back(Point(w(0),w(1)));
		pv.push_back(Point(w(2),w(3)));
		float cex=pv[0].x+x*(pv[1].x-pv[0].x)+y*(pv[0].y-pv[1].y);
		float cey=pv[0].y+x*(pv[1].y-pv[0].y)+y*(pv[1].x-pv[0].x);
		pv.push_back(Point(cex,cey));
		shadow[t]=pv;
	}

	
	VectorXf f(Points.size()*2);
	f.setZero();
	for(auto t:Triangles){
		float dn=(shadow[t][1].x-shadow[t][0].x)*(shadow[t][1].x-shadow[t][0].x)+(shadow[t][1].y-shadow[t][0].y)*(shadow[t][1].y-shadow[t][0].y);
		float dr=(t->OrdinaryPoints[1].x-t->OrdinaryPoints[0].x)*(t->OrdinaryPoints[1].x-t->OrdinaryPoints[0].x)+(t->OrdinaryPoints[1].y-t->OrdinaryPoints[0].y)*(t->OrdinaryPoints[1].y-t->OrdinaryPoints[0].y);
		float r=sqrt(dr/dn);
//		cout<<r<<endl;
		for(int p=0;p<3;p++){
			int i=PointsOrder[t->v[(p+0)%3]];
			int j=PointsOrder[t->v[(p+1)%3]];
			float ex=r*(shadow[t][(p+1)%3].x-shadow[t][(p+0)%3].x);
			float ey=r*(shadow[t][(p+1)%3].y-shadow[t][(p+0)%3].y);
			f(i*2+0)+=2*ex;
			f(i*2+1)+=2*ey;
			f(j*2+0)+=-2*ex;
			f(j*2+1)+=-2*ey;
		}
	}
	MatrixXf HM01(un*2,qn*2);
	MatrixXf HM10(qn*2,un*2);
	HM01=EdgesTopologyM.block(0, un*2, un*2, qn*2);
	HM10=EdgesTopologyM.block(un*2, 0, qn*2, un*2);
	MatrixXf D=HM01+HM10.transpose();
//	cout<<endl<<FinallyM<<endl<<D<<endl;
	VectorXf f0=f.block(0, 0, un*2, 1);
	u.setZero();
	u=-FinallyM*(D*q+f0);
	
	for(int i=0;i<un;i++){
		Points[i]->x=u(i*2+0);
		Points[i]->y=u(i*2+1);
	}
}






