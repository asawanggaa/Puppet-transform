//
//  main.cpp
//  Puppet-transform
//
//  Created by t_wangju on 03/02/2017.
//  Copyright © 2017 t_wangju. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <GLUT/GLUT.h>
#include <OpenGL/gl.h>

#include <opencv2/opencv.hpp>
#include <eigen3/Eigen/Dense>

#include "shape.hpp"
#include "TestFunc.hpp"
#include "InitFunc.hpp"
#include "ptran.hpp"

using namespace std;
using namespace Eigen;
using namespace psl2t;

float theta(Point* p,Point *i,Point *j){
	float xi=i->x-p->x;
	float yi=i->y-p->y;
	float xj=j->x-p->x;
	float yj=j->y-p->y;
	float dot=xi*xj+yi*yj;
	float cross=xi*yj-xj*yi;
	return cvFastArctan(cross, dot);
}

bool  InPoly(vector<Point*> input, Point* p){
	float th=0;
	for(auto i=input.begin();i!=input.end();i++){
		auto j=i+1;
		if(j==input.end()) j=input.begin();
		float t=theta(p, *i, *j);
		if(t>180) t=t-360;
//		cout<<t<<endl;
		th+=t;
	}
	if(th>-180&&th<180) return true;
	return false;
}

GLuint tex;
string texpath("/Users/t_wangju/Workplace/chessboard.jpg");
string datapath("/Users/t_wangju/Workplace/Puppet-transform/dude.dat");

ptran *tp;
set<Triangle*> res;

void DisplayFunc(void){
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.5, 0.5, 0.5);
	glEnable(GL_TEXTURE_2D);

	for(auto t:res){
		glBegin(GL_LINE_LOOP);
		glVertex2f(t->v[0]->x, t->v[0]->y);
		glVertex2f(t->v[1]->x, t->v[1]->y);
		glVertex2f(t->v[2]->x, t->v[2]->y);
		glEnd();
	}
	
	glColor3f(1.0f, 0.0f, 0.0f);
	for(auto p:tp->ConstraintPoints){
		DrawHandle(p->x, p->y, 5);
	}

	glFlush();
	glutSwapBuffers();
}

Point *poi=NULL;;

void MouseClickFunc(int button, int state, int x, int y){
	if(button==GLUT_LEFT_BUTTON&&state==GLUT_DOWN){
//		glClear(GL_COLOR_BUFFER_BIT);
		y=800-y;
		
		for(auto p:tp->ConstraintPoints){
			float tx=x-p->x;
			float ty=y-p->y;
			if(tx*tx+ty*ty<25){
				poi=p;
			}
		}
		
//		glColor3f(0.0f, 0.0f, 0.0f);
//		for(auto t:res){
//			glBegin(GL_LINE_LOOP);
//			glVertex2f(t->v[0]->x, t->v[0]->y);
//			glVertex2f(t->v[1]->x, t->v[1]->y);
//			glVertex2f(t->v[2]->x, t->v[2]->y);
//			glEnd();
//		}
		
//		glColor3f(1.0f, 0.0f, 0.0f);
//		for(auto p:tp->ConstraintPoints){
//			DrawHandle(p->x, p->y, 5);
//		}
//		glFlush();
//		glutSwapBuffers();
	}
	if(button==GLUT_LEFT_BUTTON&&state==GLUT_UP){
		y=800-y;
		if(!poi) return;
		poi->x=x;
		poi->y=y;
//		tp->FlushConstriant();
//		glColor3f(0.0f, 0.0f, 0.0f);
//		for(auto t:res){
//			glBegin(GL_LINE_LOOP);
//			glVertex2f(t->v[0]->x, t->v[0]->y);
//			glVertex2f(t->v[1]->x, t->v[1]->y);
//			glVertex2f(t->v[2]->x, t->v[2]->y);
//			glEnd();
//		}
//
//		glColor3f(1.0f, 0.0f, 0.0f);
//		for(auto p:tp->ConstraintPoints){
//			DrawHandle(p->x, p->y, 5);
//		}
//		poi=NULL;
//		glFlush();
//		glutSwapBuffers();
		poi=NULL;
	}
}
void MouseMotionFunc(int x, int y){
	if(!poi) return;
	glClear(GL_COLOR_BUFFER_BIT);
	y=800-y;
	poi->x=x;
	poi->y=y;
	tp->FlushConstriant();
	glColor3f(0.0f, 0.0f, 0.0f);
	for(auto t:res){
		glBegin(GL_LINE_LOOP);

		glVertex2f(t->v[0]->x, t->v[0]->y);
		glVertex2f(t->v[1]->x, t->v[1]->y);
		glVertex2f(t->v[2]->x, t->v[2]->y);
		glEnd();
	}
	
	glColor3f(1.0f, 0.0f, 0.0f);
	for(auto p:tp->ConstraintPoints){
		DrawHandle(p->x, p->y, 5);
	}
	glFlush();
	glutSwapBuffers();
}
int main(int argc, char * argv[]) {
	// insert code here...
	//	datainit(10, 20);
	init(argc, argv, 0, 0, 800, 800);
	tex=texinit(texpath.c_str());
	
	vector<Point*> input;
	ifstream in(datapath);
	string line;
	while(getline(in, line)){
		float a,b;
		sscanf(line.c_str(),"%f %f",&a,&b);
		input.push_back(new Point(a,b));
	}

//	input.clear();
//	input.push_back(new Point(200,400));
//	input.push_back(new Point(250,250));
//	input.push_back(new Point(400,200));
//	input.push_back(new Point(550,250));
//	input.push_back(new Point(600,400));
//	input.push_back(new Point(550,550));
//	input.push_back(new Point(400,600));
//	input.push_back(new Point(250,550));
	
	Mesh m(100, 100, 700, 700);
	auto t1=new Point(350,470);
	auto t2=new Point(330,500);
	auto t3=new Point(360,380);
	auto tinput=input;
	tinput.push_back(t1);
	tinput.push_back(t2);
	tinput.push_back(t3);
	
	res=m.Bowyer_Watson(tinput);
	res=m.EraseEdgeTriangles();
	
	
	set<Triangle*> bts;
	for(auto t:res){
		float x=0,y=0;
		for(auto p:t->v){
			x+=p->x;
			y+=p->y;
		}
		Point mid(x/3,y/3);
		if(InPoly(input, &mid)){
			bts.insert(t);
		}
	}
	for(auto t:bts){
		res.erase(t);
		delete t;
	}
	
//	tp=new puppet(res);
//	tp->InsertConstraint(t1);
//	tp->InsertConstraint(t2);
//	tp->InsertConstraint(t3);
	tp=new ptran(res);
	
//	for(int i=1;i<tinput.size();i++){
//		if(i%20==0){
//			tp->InsertConstraintPoints(tinput[i]);
//		}
//	}
	
	tp->InsertConstraintPoints(t1);
	tp->InsertConstraintPoints(t2);
	tp->InsertConstraintPoints(t3);
	tp->SetConstriant();
	tp->FlushConstriant();
//	tp->flush();
//
//	MatrixXf G=tp->GMatrix.block(0, 0, 8, 8);
//	cout<<endl<<G<<endl;
//	MatrixXf H;
//	H=G.transpose();
//	G=G+H;
//	cout<<endl<<G<<endl;
//	cout<<endl<<G.inverse()<<endl;
	
	glutDisplayFunc(DisplayFunc);
	
	glutMouseFunc(MouseClickFunc);
	glutMotionFunc(MouseMotionFunc);
//	glutKeyboardFunc(KeyboardFunc);
	
	glutMainLoop();
	return 0;
}
