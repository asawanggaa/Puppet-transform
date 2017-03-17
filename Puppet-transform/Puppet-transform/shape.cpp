//
//  shape.cpp
//  Puppet-transform
//
//  Created by t_wangju on 03/02/2017.
//  Copyright © 2017 t_wangju. All rights reserved.
//

#include "shape.hpp"

using namespace std;

namespace psl2t {
	Point::Point(float ix,float iy):x(ix),y(iy){}
	Point::~Point(){}
	int Point::left(Point *A, Point *B)const{
		Eigen::Matrix2f M;
		M<<A->x-x,A->y-y,B->x-x,B->y-y;
		float ans=M.determinant();
		if(ans>1e-8){
			return 1;
		}else if(ans<-1e-8){
			return -1;
		}else{
			return 0;
		}
	}
	float Point::dis(const Point &A)const{
		return sqrt((x-A.x)*(x-A.x)+(y-A.y)*(y-A.y));
	}
	
	Triangle::Triangle(Point *l,Point *m,Point *n):O(0,0){
		v[0]=l,v[1]=m,v[2]=n;
		if(v[2]->left(v[0], v[1])==-1) swap(v[0], v[1]);
		for(auto p:v) p->linkT.insert(this);
		Eigen::Matrix3f A,B,C;
		A<<
		l->x*l->x+l->y*l->y,l->y,1,
		m->x*m->x+m->y*m->y,m->y,1,
		n->x*n->x+n->y*n->y,n->y,1;
		B<<
		l->x,l->x*l->x+l->y*l->y,1,
		m->x,m->x*m->x+m->y*m->y,1,
		n->x,n->x*n->x+n->y*n->y,1;
		C<<
		l->x,l->y,1,
		m->x,m->y,1,
		n->x,n->y,1;
		O.x=A.determinant()/2/C.determinant();
		O.y=B.determinant()/2/C.determinant();
		R=0;
		for(auto p:v) R+=(p->dis(O)/3.0f);
		
		for(int i=0;i<3;i++){
			Point *s=v[(i+1)%3];
			Point *t=v[(i+2)%3];
			
			coeffi[i][0]=((v[i]->x-s->x)*(t->x-s->x)+(v[i]->y-s->y)*(t->y-s->y))/((t->x-s->x)*(t->x-s->x)+(t->y-s->y)*(t->y-s->y));
			coeffi[i][1]=((v[i]->x-s->x)*(s->y-t->y)+(v[i]->y-s->y)*(t->x-s->x))/((t->x-s->x)*(t->x-s->x)+(t->y-s->y)*(t->y-s->y));
			OrdinaryPoints.push_back(Point(v[i]->x,v[i]->y));
		}
	}
	Triangle::~Triangle(){
		for(auto p:v) p->linkT.erase(this);
	}
	bool Triangle::in(const Point *p)const{
		vector<int>	res;
		for(int i=0;i<3;i++){
			Point *A=v[(i+0)%3];
			Point *B=v[(i+1)%3];
			Point *C=v[(i+2)%3];
			res.push_back(p->left(A,B)*C->left(A,B));
		}
		for(auto i:res) if(i<-1e-8) return false;
		return true;
	}
	bool Triangle::incircle(const Point *p)const{
		if(p->dis(O)<R+1e-8){
			return true;
		}else{
			return false;
		}
	}
	Triangle* Triangle::next(const Point* p)const{
		for(int i=0;i<3;i++){
			Point *A=v[(i+0)%3];
			Point *B=v[(i+1)%3];
			Point *C=v[(i+2)%3];
			if(p->left(A,B)*C->left(A,B)<0&&p->left(A,C)*B->left(A,C)>=0){
				for(auto t:A->linkT){
					if(B->linkT.count(t)){
						if(t!=this){
							return t;
						}
					}
				}
			}
		}
		cerr<<"error:not found a neibour edge;"<<endl;
		exit(1);
	}
	
	set<Triangle*> Mesh::InsertAPoint(set<Triangle*> &ts, Point *p){
		Triangle *fbt=FindFirstBadTriangle(ts, p);
		set<Triangle*> badts=ExploreBadTriangles(ts, fbt, p);
		set<set<Point*>> use=FindUnSharedEdges(badts);
		for(auto t:badts) ts.erase(t),delete t;
		set<Triangle*> nts=ConstructNewTriangles(use, p);
		for(auto t:nts) ts.insert(t);
		return ts;
	}
	Triangle* Mesh::FindFirstBadTriangle(set<Triangle*> &ts,Point *p){
		Triangle* t=*ts.begin();
		while(!t->incircle(p)){
			t=t->next(p);
		}
		return t;
	}
	set<Triangle*> Mesh::ExploreBadTriangles(set<Triangle*> &ts,Triangle *fbt,Point *p){
		set<Triangle*> badts;
		queue<Triangle*> wait;
		set<Triangle*> visited;
		wait.push(fbt);
		while(!wait.empty()){
			Triangle *t=wait.front();
			wait.pop();
			if(visited.count(t)){
				continue;
			}
			visited.insert(t);
			if(!t->incircle(p)){
				continue;
			}
			badts.insert(t);
			for(int i=0;i<3;i++){
				Point* A=t->v[(i+0)%3];
				Point* B=t->v[(i+1)%3];
				for(auto st:A->linkT){
					if(B->linkT.count(st)){
						if(st!=t){
							wait.push(st);
						}
					}
				}
			}
		}
		return badts;
	}
	set<set<Point*>> Mesh::FindUnSharedEdges(set<Triangle *> &badts){
		set<set<Point*>> res;
		for(auto t:badts){
			for(int i=0;i<3;i++){
				Point* A=t->v[(i+0)%3];
				Point* B=t->v[(i+1)%3];
				set<Point*> e;
				e.insert(A);
				e.insert(B);
				if(res.count(e)){
					res.erase(e);
				}else{
					res.insert(e);
				}
			}
		}
		return res;
	}
	set<Triangle*> Mesh::ConstructNewTriangles(set<set<Point *>> &use,Point *p){
		set<Triangle*> res;
		for(auto e:use){
			if(e.size()!=2){
				cerr<<"an edge has more than 2 points"<<endl;
				exit(2);
			}
			Triangle *t=new Triangle(*e.begin(),*e.rbegin(),p);
			res.insert(t);
		}
		return res;
	}
	set<Triangle*> Mesh::Bowyer_Watson(vector<Point *> &input){
		for(auto p:input){
			ts=InsertAPoint(ts, p);
//			cout<<p->x<<' '<<p->y<<endl;
		}
		return ts;
	}
	set<Triangle*> Mesh::EraseEdgeTriangles(){
		set<Triangle*> rubbish;
		for(auto t:LB->linkT) rubbish.insert(t);
		for(auto t:RB->linkT) rubbish.insert(t);
		for(auto t:LT->linkT) rubbish.insert(t);
		for(auto t:RT->linkT) rubbish.insert(t);
		for(auto t:rubbish){
			ts.erase(t);
			delete t;
		}
		return ts;
	}

	set<Triangle*> Mesh::FindCrossTriangles(set<Triangle*> &ts,Point* s,Point* t){
		set<Triangle*> res;
		if(s==t) return res;
		
		Point *A=NULL,*B=NULL,*C=NULL;
		Triangle *st=NULL;
		// find out st & a,b;
		for(auto tt:s->linkT){
			for(int i=0;i<3;i++){
				A=tt->v[(0+i)%3];
				B=tt->v[(1+i)%3];
				if(A!=s&&B!=s)
					break;
			}
			if(A->left(s, t)*B->left(s, t)<0&&s->left(A, B)*t->left(A, B)<0){
				st=tt;
				if(A->left(s, t)<0) swap(A, B);
				break;
			}
		}
		
		res.insert(st);
		while(C!=t){
			for(auto tt:A->linkT){
				if(B->linkT.count(tt)&&tt!=st){
					st=tt;
					break;
				}
			}
			for(auto p:st->v) if(p!=A&&p!=B) C=p;
			res.insert(st);
			switch (C->left(s, t)) {
				case -1:
					B=C;
					break;
				case 1:
					A=C;
					break;
				case 0:
					cerr<<"A point is on this segment"<<endl;
					break;
			}
			
		}
		return res;
	}
	
	Mesh::Mesh(float xmin,float ymin,float xmax,float ymax){
		LB=new Point(xmin,ymin);
		LT=new Point(xmin,ymax);
		RB=new Point(xmax,ymin);
		RT=new Point(xmax,ymax);
		ts.insert(new Triangle(LB,RT,LT));
		ts.insert(new Triangle(LB,RT,RB));
	}
	
	
}
