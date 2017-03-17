//
//  InitFunc.hpp
//  Algorithm Test
//
//  Created by t_wangju on 23/01/2017.
//  Copyright © 2017 t_wangju. All rights reserved.
//

#ifndef InitFunc_hpp
#define InitFunc_hpp

#include <iostream>
#include <GLUT/GLUT.h>
#include <OpenGL/OpenGL.h>
#include <opencv2/opencv.hpp>

void init(int argc, char * argv[], int x, int y, int width, int height);
GLuint texinit(const char * path);

void init(int argc, char * argv[], int x, int y, int width, int height){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA);
	glutInitWindowPosition(x, y);
	glutInitWindowSize(width, height);
	glutCreateWindow("Puppet Transform");
	gluOrtho2D(x, x+width, y, y+height);
	
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glReadBuffer(GL_FRONT);
	glDrawBuffer(GL_BACK);
	glFlush();
	glutSwapBuffers();
}
GLuint texinit(const char * path){
	GLuint tex;
	glShadeModel(GL_SMOOTH);
	glEnable(GL_COLOR_MATERIAL);
	glMatrixMode(GL_PROJECTION);
	
	cv::Mat img = cv::imread(path);
	cv::imshow("chessboard", img);
	cv::flip(img, img, 0);
	
	glPixelStorei(GL_UNPACK_ALIGNMENT, (img.step&3)?1:4);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, GLint(img.step/img.elemSize()));
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_2D, tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img.cols, img.rows, 0, GL_BGR, GL_UNSIGNED_BYTE, img.ptr());
	glGenerateMipmap(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);
	return tex;
}

#endif /* InitFunc_hpp */
