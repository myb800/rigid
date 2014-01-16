// rigid_body.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "rigid.h"
int wind_h=600,wind_w=600;
squares *sq;
number_t mouse_x=0,mouse_y=0;

void display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0,1.0,1.0);
	vector<rigid_sys::body_info> toshow=sq->current_body();
	
	for(int i=0;i<(int)toshow.size();i++)
	{
		glColor3f(0.0,i,0);
		glBegin(GL_LINE_LOOP);
			for(int j=0;j<4;j++)
				glVertex3f(toshow[i].pos[j].x+sq->dynamics[i].x.x,toshow[i].pos[j].y+sq->dynamics[i].x.y,toshow[i].pos[j].z);
		glEnd();
	}
	//arrow
	number_t pi=3.1415926;
	number_t theta=15*pi/180;
	number_t length=0.1;
	glBegin(GL_LINES);
		glVertex3f(-1,-1,0);
		glVertex3f(mouse_x,mouse_y,0);
	glEnd();

	number_t mouse_angle=pi+atan((mouse_y+1)/(mouse_x+1));

	glBegin(GL_LINES);
		glVertex3f(mouse_x+length*cos(mouse_angle+theta),mouse_y+length*sin(mouse_angle+theta),0);
		glVertex3f(mouse_x,mouse_y,0);
		glVertex3f(mouse_x+length*cos(mouse_angle-theta),mouse_y+length*sin(mouse_angle-theta),0);
		glVertex3f(mouse_x,mouse_y,0);
	glEnd();

	glutSwapBuffers();
}
void init()
{
	glClearColor(1.0,1.0,1.0,1.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1.0,1.0,-1.0,1.0,-1.0,1.0);

	vector<rigid_sys::body_info> vbi;
	vector<rigid_sys::dynamic_info> vdi;
	rigid_sys::body_info bi;
	bi.total_mass=4;
	bi.pos.push_back(vec(0.1,0.1,0))  ;bi.pos.push_back(vec(0.1,-0.1,0));
	bi.pos.push_back(vec(-0.1,-0.1,0));bi.pos.push_back(vec(-0.1,0.1,0));
	bi.mass.push_back(1);bi.mass.push_back(1);
	bi.mass.push_back(1);bi.mass.push_back(1);
	bi.generate_Ibody_inv();
	rigid_sys::dynamic_info di;
	di.L=vec(0,0,0);di.P=vec(0,0,0);di.q=quaternion(1,0,0,0);di.x=vec(0.5,0,0);
	vbi.push_back(bi);vdi.push_back(di);

	bi.total_mass=4;bi.pos.clear();bi.mass.clear();
	bi.pos.push_back(vec(0.1,0.1,0))  ;bi.pos.push_back(vec(0.1,-0.1,0));
	bi.pos.push_back(vec(-0.1,-0.1,0));bi.pos.push_back(vec(-0.1,0.1,0));
	bi.mass.push_back(1);bi.mass.push_back(1);
	bi.mass.push_back(1);bi.mass.push_back(1);
	bi.generate_Ibody_inv();
	di.L=vec(0,0,0.01);di.P=vec(0,0,0);di.q=quaternion(0.8,0,0,0.6);di.x=vec(-0.5,0,0);
	vbi.push_back(bi);vdi.push_back(di);

	sq=new squares(vbi,vdi);
	
}
void timerf(int value)
{
	number_t h=1;
	sq->update(value,value+h);
	glutPostRedisplay();
	glutTimerFunc(20,timerf,value+h);
}
void reshape(int w,int h)
{
	wind_h=h;
	wind_w=w;
    glColor3f(1.0,1.0,1.0);
    glViewport(0,0,(GLsizei)w,(GLsizei)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	gluOrtho2D(-1.0,1.0,-1.0,1.0);
    glutSwapBuffers();
}
void motion(int x,int y)
{
	mouse_x=2*(number_t)x/wind_w-1;
	mouse_y=1-2*(number_t)y/wind_h;
}
void mouse(int b,int s,int x,int y)
{
	if(b==GLUT_LEFT_BUTTON && s==GLUT_DOWN)
	{
		rigid_sys::body_info bi;
		bi.total_mass=4;bi.pos.clear();bi.mass.clear();
		bi.pos.push_back(vec(0.1,0.1,0))  ;bi.pos.push_back(vec(0.1,-0.1,0));
		bi.pos.push_back(vec(-0.1,-0.1,0));bi.pos.push_back(vec(-0.1,0.1,0));
		bi.mass.push_back(1);bi.mass.push_back(1);
		bi.mass.push_back(1);bi.mass.push_back(1);
		bi.generate_Ibody_inv();
		sq->square_flyin(bi,vec((1+mouse_x)/50,(1+mouse_y)/50,0),vec(-1,-1,0),1);

	}
}
int _tmain(int argc, _TCHAR* argv[])
{

	glutInit(&argc,(char**)argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(wind_w,wind_h);
	glutInitWindowPosition(250,250);
	glutCreateWindow("Square Collision");

	init();
	glutDisplayFunc(&display);
	glutReshapeFunc(&reshape);
	glutMouseFunc(&mouse);
	glutTimerFunc(20,timerf,0);
	//glutKeyboardFunc(&keyboard);
	glutPassiveMotionFunc(&motion);
	glutMainLoop();
	getchar();
	return 0;
}

