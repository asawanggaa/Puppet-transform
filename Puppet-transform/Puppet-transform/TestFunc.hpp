//
//  TestFunction.hpp
//  LinearWarp
//
//  Created by t_wangju on 03/02/2017.
//  Copyright © 2017 t_wangju. All rights reserved.
//

#ifndef TestFunction_hpp
#define TestFunction_hpp

#include <iostream>
#include <GLUT/GLUT.h>
#include <OpenGL/OpenGL.h>
#include <opencv2/opencv.hpp>

#define HandleRadius 5

void DrawTex(GLuint tex, float x, float y, float width, float height);
void DrawHandle(int x, int y, int Radius){
	double theta;
	
	glBegin(GL_LINE_LOOP);
	{
		for (int i=0; i<36; i++) {
			theta=i*2*M_PI/36;
			double dx=Radius*cos(theta)+x;
			double dy=Radius*sin(theta)+y;
			glVertex2d(dx, dy);
		}
	}
	glEnd();
}
void DrawTex(GLuint tex, float x, float y, float width, float height){
	glDisable(GL_BLEND);
	glReadBuffer(GL_FRONT);
	glDrawBuffer(GL_BACK);
	
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1, 1, 1, 1);
	glColor3f(1.0, 1.0, 1.0);
	glBindTexture(GL_TEXTURE_2D, tex);
	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
	{
		glTexCoord2f(0.0f, 0.0f);
		glVertex2f(x, y);
		glTexCoord2f(1.0f, 0.0f);
		glVertex2f(x+width, y);
		glTexCoord2f(1.0f, 1.0f);
		glVertex2f(x+width, y+height);
		glTexCoord2f(0.0f, 1.0f);
		glVertex2f(x, y+height);
	}
	glEnd();
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);
	
}
#endif /* TestFunction_hpp */
