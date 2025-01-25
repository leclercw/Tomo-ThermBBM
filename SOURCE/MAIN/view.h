#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <time.h> 
#include <sys/time.h> 
#include <sys/resource.h> 
#include <string.h>
#include <GL/gl.h>
#include <GL/glu.h> /* fichiers d?entetes OpenGL, GLU et GLUT */
#include <GL/freeglut.h>
#include <map>
#include <cassert>
#include <vector>
#include <limits>

//#include "view.h"

typedef double R;

using namespace std;


////////////////////////////////////////////////////////////////////////////////
//// Paramètres du système
int cH,nQ,H1,H2,H3,HX,HY,HZ,NB_SPH,NB_DIS,NB_FIB,N_TOT,V_TOT;
R pC,cont,PR_INC,R_DIS,R_ECH,FACF;

R minsxx,minsyy,minszz,minsxy,minsxz,minsyz;
R maxsxx,maxsyy,maxszz,maxsxy,maxsxz,maxsyz;
R mintt,minux,minuy,minuz;
R maxtt,maxux,maxuy,maxuz;
R minfx,minfy,minfz;
R maxfx,maxfy,maxfz;

////////////////////////////////////////////////////////////////////////////////
//// Tables 

bool * LIST_N;
int  * LIST_A;

R * LIST_PS;
R * LIST_GA;
R * LIST_PH;

int ** LIST_Q;
int * LIST_V;

R    * LIST_R;
R    * LIST_X;
R    * LIST_Y;
R    * LIST_Z;

R    * LIST_SXX;
R    * LIST_SXY;
R    * LIST_SXZ;
R    * LIST_SYY;
R    * LIST_SYZ;
R    * LIST_SZZ;
R    * LIST_UX;
R    * LIST_UY;
R    * LIST_UZ;

R    * LIST_T;
R    * LIST_FX;
R    * LIST_FY;
R    * LIST_FZ;

////////////////////////////////////////////////////////////////////////////////  
//// Paramètres visu

float theta_eye;
float phi_eye;
float rayon;
float x_eye;
float y_eye;
float z_eye;
int x_g;
int x_d;
int y_g;
int y_d;
int click_g;
int click_d;

int cF=1;
int posi,dire;
int nb;
bool running[1000000];

////////////////////////////////////////////////////////////////////////////////  


void InitGL2(void)
{
glClearColor(0,0,0,0);
/*definition de la couleur utilisee */
/*pour effacer */
glColor3f(1.0,1.0,1.0);

nb=0;
}

void InitGL3(void)
{
glClearColor(0,0,0,0);
/*definition de la couleur utilisee */
/*pour effacer */
glColor3f(1.0,1.0,1.0);

	// remise à 0 si space
	theta_eye=0;
	phi_eye=0;	
	rayon=2.5*((GLdouble) cH)/((GLdouble) cF);
	x_eye=0;
	y_eye=0;
	z_eye=rayon;
	x_g=280;
	y_g=280;
	x_d=280;
	y_d=280;

nb=0;
dire=1;
posi=cH/2;
}


void setup_illumination(){
  // Intialise and set lighting parameters
  GLfloat light_pos[] = {0.,0.,0., 1.};
  GLfloat light_ka[] = {0., 0., 0.2, 1.};
  GLfloat light_kd[] = {1., 1., 1., 1.};
  GLfloat light_ks[] = {1., 1., 1., 1.};

  glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
  glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ka);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_kd);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_ks);

  
  // Initialise and set material parameters
  GLfloat material_ka[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat material_kd[] = {0.43, 0.47, 0.54, 1.0};
  GLfloat material_ks[] = {0.33, 0.33, 0.52, 1.0};
  GLfloat material_ke[] = {0.0, 0.0, 0.0, 0.0};
  GLfloat material_se[] = {10.0};

  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  material_ka);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  material_kd);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  material_ks);
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION,  material_ke);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, material_se);
  
// Smooth shading
glShadeModel(GL_SMOOTH);

//Enable lighting
glEnable (GL_LIGHTING);
glEnable (GL_LIGHT0);

}

void afficherColorbar_t()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = mintt;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxtt;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficherColorbar_phix()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minfx;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxfx;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}


void afficherColorbar_phiy()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minfy;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxfy;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}


void afficherColorbar_phiz()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minfz;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxfz;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficherColorbar_sxx()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minsxx;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxsxx;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficherColorbar_sxy()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minsxy;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxsxy;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficherColorbar_sxz()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minsxz;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxsxz;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficherColorbar_syy()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minsyy;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxsyy;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficherColorbar_syz()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minsyz;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxsyz;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficherColorbar_szz()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minszz;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxszz;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficherColorbar_ux()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150.,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minux;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8.,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxux;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8.,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}


void afficherColorbar_uy()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minuy;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxuy;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficherColorbar_uz()
{

glClearDepth(1.);
glClearColor(1.0,1.0,1.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(0,0,150,0.,0.,0.,0,1,0);

//Colorbar values
glColor3f(0.0f,0.0f,0.0f);
float fVal ;
fVal = minuz;
char cVal[32] ;
sprintf(cVal,"%.2E",fVal) ;

cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,-50.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

glColor3f(0.0f,0.0f,0.0f);
fVal = maxuz;
sprintf(cVal,"%.2E",fVal) ;
cVal[9]=' '  ;
cVal[10]='P'  ;
cVal[11]='a'  ;
glRasterPos2f(-8,50.-4.);

for (int k=0;k<=11;k++ ) {
glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,cVal[k]);
}

//colorbar scale
   
    for(int i=0;i<100;i++){
		
	glColor3f(1.,(1.-i/100.),0.); 	
	glBegin(GL_QUADS);
	glVertex2f(-20,i-50);
	glVertex2f(-10,i-50);
	glVertex2f(-10,(i+1)-50);
	glVertex2f(-20,(i+1)-50);	
    glEnd(); 
    }

// cadre colorbar

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f(-20,50);
	glVertex2f(-10,50);
	glVertex2f(-10,-50);
	glVertex2f(-20,-50);		  
   	glEnd();

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher2_agr()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Microstructure";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

//infos
/*glColor3f(255.,255.,255.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string stringnC = static_cast<ostringstream*>( &(ostringstream() << nC) )->str();
string sstrnF = "Nombre d inclusions : " + stringnC + "";
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

glColor3f(255.,255.,255.); 
glRasterPos2f(-0.5*cH/cF,-0.6*cH/cF);
string stringpC0 = static_cast<ostringstream*>( &(ostringstream() << pC0) )->str();
string stringpC = static_cast<ostringstream*>( &(ostringstream() << pC) )->str();
string sstrpF = "Pourcentage d inclusions : " + stringpC + "/"+ stringpC + "" ;
const char * cstrpF = sstrpF.c_str();
const unsigned char *info2 = reinterpret_cast< const unsigned char * >(cstrpF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info2 ); */
long nbe = 0;

R col1[7];
R col2[7];
R col3[7];

col1[0]=1;col1[1]=0;col1[2]=1;col1[3]=1;col1[4]=0;col1[5]=0;col1[6]=0.5;
col2[0]=1;col2[1]=1;col2[2]=0;col2[3]=0;col2[4]=1;col2[5]=0;col2[6]=0.5;
col3[0]=0;col3[1]=1;col3[2]=1;col3[3]=0;col3[4]=0;col3[5]=1;col3[6]=0.5;


for(int j=0;j<cH;j++){ 
  for(int i=0;i<cH;i++){

    glColor3f((LIST_A[nbe]>0)?(col1[LIST_A[nbe]%7]):1,(LIST_A[nbe]>0)?(col2[LIST_A[nbe]%7]):1,(LIST_A[nbe]>0)?(col3[LIST_A[nbe]%7]):1);  
    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 /*   
    glColor3f(0.,0.,0.); 
    glBegin(GL_LINE_LOOP);
    
	glVertex2i((i-cH/2)/cF,(j-cH/2)/cF);
	glVertex2i((i+1-cH/2)/cF,(j-cH/2)/cF);
	glVertex2i((i+1-cH/2)/cF,(j+1-cH/2)/cF);
	glVertex2i((i-cH/2)/cF,(j+1-cH/2)/cF);
	
    glEnd(); */     
   nbe++;
 }  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}


void afficher2_cont()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Microstructure";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

//infos
/*glColor3f(255.,255.,255.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string stringnC = static_cast<ostringstream*>( &(ostringstream() << nC) )->str();
string sstrnF = "Nombre d inclusions : " + stringnC + "";
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

glColor3f(255.,255.,255.); 
glRasterPos2f(-0.5*cH/cF,-0.6*cH/cF);
string stringpC0 = static_cast<ostringstream*>( &(ostringstream() << pC0) )->str();
string stringpC = static_cast<ostringstream*>( &(ostringstream() << pC) )->str();
string sstrpF = "Pourcentage d inclusions : " + stringpC + "/"+ stringpC + "" ;
const char * cstrpF = sstrpF.c_str();
const unsigned char *info2 = reinterpret_cast< const unsigned char * >(cstrpF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info2 ); */
long nbe = 0;

for(int j=0;j<cH;j++){ 
  for(int i=0;i<cH;i++){

    glColor3f(255.,(1-LIST_N[nbe])*255.,0.);  
    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 /*   
    glColor3f(0.,0.,0.); 
    glBegin(GL_LINE_LOOP);
    
	glVertex2i((i-cH/2)/cF,(j-cH/2)/cF);
	glVertex2i((i+1-cH/2)/cF,(j-cH/2)/cF);
	glVertex2i((i+1-cH/2)/cF,(j+1-cH/2)/cF);
	glVertex2i((i-cH/2)/cF,(j+1-cH/2)/cF);
	
    glEnd(); */     
   nbe++;
 }  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher2_disk()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Microstructure";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

//infos
/*glColor3f(255.,255.,255.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string stringnC = static_cast<ostringstream*>( &(ostringstream() << nC) )->str();
string sstrnF = "Nombre d inclusions : " + stringnC + "";
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

glColor3f(255.,255.,255.); 
glRasterPos2f(-0.5*cH/cF,-0.6*cH/cF);
string stringpC0 = static_cast<ostringstream*>( &(ostringstream() << pC0) )->str();
string stringpC = static_cast<ostringstream*>( &(ostringstream() << pC) )->str();
string sstrpF = "Pourcentage d inclusions : " + stringpC + "/"+ stringpC + "" ;
const char * cstrpF = sstrpF.c_str();
const unsigned char *info2 = reinterpret_cast< const unsigned char * >(cstrpF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info2 ); */

glColor3f(0.,0.,0.);
for(int k=0;k<NB_DIS;k++)
{

   glBegin(GL_LINE_STRIP);
    for(int j=0; j <= 32; j++)
    {
      
        double lng1 = double(j)*(2.0*3.14159)/double(32);
        
	double x1 = (LIST_X[k]+cos(lng1)*LIST_R[k]-cH/2)/cF;	
        double y1 = (LIST_Y[k]+sin(lng1)*LIST_R[k]-cH/2)/cF;
        glVertex2f(x1,y1);
	
	
    }
    glEnd();

}


// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 


glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}



void afficher2_t()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Champ de temperature";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

R imins2=mintt;
R imaxs2=maxtt;

long nbe = 0;

for(int j=0;j<cH;j++){ 
  for(int i=0;i<cH;i++){

    int alpha = nbe/cH;
    int beta  = nbe%cH;
    
    int nump  = alpha*cH+beta;
    R uu = LIST_T[nump];
    
    glColor3f(1.,1.-pow(((uu-mintt)/(maxtt-mintt)),1),0.); 
    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
   nbe++;
 }  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 
      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher2_phix()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Densite de flux x";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

R imins2=minfx;
R imaxs2=maxfx;

long nbe = 0;

for(int j=0;j<cH;j++){ 
  for(int i=0;i<cH;i++){

    int alpha = nbe/cH;
    int beta  = nbe%cH;
    
    int nump  = alpha*cH+beta;
    R uu = LIST_FX[nump];
    
    glColor3f(1.,1.-pow(((uu-minfx)/(maxfx-minfx)),1),0.); 
    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 

   nbe++;
 }  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 
      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}


void afficher2_phiy()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Densite de flux y";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

R imins2=minfy;
R imaxs2=maxfy;

long nbe = 0;

for(int j=0;j<cH;j++){ 
  for(int i=0;i<cH;i++){

    int alpha = nbe/cH;
    int beta  = nbe%cH;
    
    int nump  = alpha*cH+beta;
    R uu = LIST_FY[nump];
    
    glColor3f(1.,1.-pow(((uu-minfy)/(maxfy-minfy)),1),0.); 
    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
    nbe++;
 }  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 
      	

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}



void afficher2_ux()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Deplacement - dir. x";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

R imins2=minux;
R imaxs2=maxux;

long nbe = 0;

for(int j=0;j<cH;j++){ 
  for(int i=0;i<cH;i++){

    int alpha = nbe/cH;
    int beta  = nbe%cH;
    
    int nump  = alpha*cH+beta;
    R uu = LIST_UX[nump];
    
    glColor3f(1.,1.-pow(((uu-minux)/(maxux-minux)),1),0.); 
    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
   nbe++;
 }  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 
      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}






void afficher2_uy()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Deplacement - dir. y";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

R imins2=minuy;
R imaxs2=maxuy;

long nbe = 0;

for(int j=0;j<cH;j++){ 
  for(int i=0;i<cH;i++){

    int alpha = nbe/cH;
    int beta  = nbe%cH;
    
    int nump  = alpha*cH+beta;
    R uu = LIST_UY[nump];
    
    glColor3f(1.,1.-pow(((uu-minuy)/(maxuy-minuy)),1),0.); 
    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
   nbe++;
 }  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 
      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}




void afficher2_sxx()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Champ de contraintes sxx";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

R imins2=minsxx;
R imaxs2=maxsxx;

long nbe = 0;

for(int j=0;j<cH;j++){ 
  for(int i=0;i<cH;i++){

    int alpha = nbe/cH;
    int beta  = nbe%cH;
    
    int nump  = alpha*cH+beta;
    R uu = LIST_SXX[nump];
    
    glColor3f(1.,1.-pow(((uu-minsxx)/(maxsxx-minsxx)),1),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
   nbe++;
 }  
}


// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher2_sxy()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Champ de contraintes sxy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

R imins2=minsxy;
R imaxs2=maxsxy;

long nbe = 0;

for(int j=0;j<cH;j++){ 
  for(int i=0;i<cH;i++){

    int alpha = nbe/cH;
    int beta  = nbe%cH;
    
    int nump  = alpha*cH+beta;
    R uu = LIST_SXY[nump];
    
    glColor3f(1.,1.-pow(((uu-minsxy)/(maxsxy-minsxy)),1),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win

   nbe++;
 }  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 
      	

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher2_syy()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Champ de contraintes syy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

R imins2=minsyy;
R imaxs2=maxsyy;

long nbe = 0;

for(int j=0;j<cH;j++){ 
  for(int i=0;i<cH;i++){

    int alpha = nbe/cH;
    int beta  = nbe%cH;
    
    int nump  = alpha*cH+beta;
    R uu = LIST_SYY[nump];
    
    glColor3f(1.,1.-pow(((uu-minsyy)/(maxsyy-minsyy)),1),0.);  
    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win

   nbe++;
 }  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 
      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_sphere()
{
  
//glClearDepth(1.);
//glClearColor(0.0,0.0,0.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(x_eye,y_eye,z_eye,0,0,0,0,1,0);

int i,j,k,l;
i=0;
j=0;
k=0;

glColor3f(1.,1.,1.);
// Box
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();
	// BACK
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j+cH-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();      
	// LEFT
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i-cH/2)/cF,  (j-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();
	// RIGHT
	glBegin(GL_LINE_LOOP);
	glVertex3i((i+cH-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF, (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();   	
	// BOTTOM
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);    
      	glEnd();

setup_illumination();
 
for(l=0;l<NB_SPH;l++){ 
// cout<<"l:"<<l<<","<<LIST_X[l]<<","<<LIST_Y[l]<<","<<LIST_Z[l]<<","<<LIST_R[l]<<endl;
  glPushMatrix();
  glTranslated((LIST_X[l]-cH/2)/cF,(LIST_Y[l]-cH/2)/cF,(LIST_Z[l]-cH/2)/cF);
  glutSolidSphere(LIST_R[l],50,50);
  glPopMatrix();
} 

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_surf()
{
  
//glClearDepth(1.);
//glClearColor(0.0,0.0,0.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice

gluLookAt(x_eye,y_eye,z_eye,0,0,0,0,1,0);

int i,j,k,l;
i=0;
j=0;
k=0;

glColor3f(1.,1.,1.);
// Box
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();
	// BACK
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j+cH-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();      
	// LEFT
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i-cH/2)/cF,  (j-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();
	// RIGHT
	glBegin(GL_LINE_LOOP);
	glVertex3i((i+cH-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF, (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();   	
	// BOTTOM
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);    
      	glEnd();
	// UP
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k+cH-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd(); 

for(l=0;l<nQ;l++){ 
 
  
	glColor3f(255.,0.,0.);  
	glBegin(GL_QUADS);
	//glNormal3i(1,0,0);
	glVertex3i(((LIST_Q[0][l]%(cH+1))-cH/2)/cF,(((LIST_Q[0][l]%((cH+1)*(cH+1)))/(cH+1))-cH/2)/cF,((LIST_Q[0][l]/((cH+1)*(cH+1)))-cH/2)/cF);  
	glVertex3i(((LIST_Q[1][l]%(cH+1))-cH/2)/cF,(((LIST_Q[1][l]%((cH+1)*(cH+1)))/(cH+1))-cH/2)/cF,((LIST_Q[1][l]/((cH+1)*(cH+1)))-cH/2)/cF);
	glVertex3i(((LIST_Q[2][l]%(cH+1))-cH/2)/cF,(((LIST_Q[2][l]%((cH+1)*(cH+1)))/(cH+1))-cH/2)/cF,((LIST_Q[2][l]/((cH+1)*(cH+1)))-cH/2)/cF);
	glVertex3i(((LIST_Q[3][l]%(cH+1))-cH/2)/cF,(((LIST_Q[3][l]%((cH+1)*(cH+1)))/(cH+1))-cH/2)/cF,((LIST_Q[3][l]/((cH+1)*(cH+1)))-cH/2)/cF);
	glEnd();
	
	glColor3f(1.,1.,1.);  
	glBegin(GL_LINE_LOOP);
	//glNormal3i(1,0,0);
	glVertex3i(((LIST_Q[0][l]%(cH+1))-cH/2)/cF,(((LIST_Q[0][l]%((cH+1)*(cH+1)))/(cH+1))-cH/2)/cF,((LIST_Q[0][l]/((cH+1)*(cH+1)))-cH/2)/cF);  
	glVertex3i(((LIST_Q[1][l]%(cH+1))-cH/2)/cF,(((LIST_Q[1][l]%((cH+1)*(cH+1)))/(cH+1))-cH/2)/cF,((LIST_Q[1][l]/((cH+1)*(cH+1)))-cH/2)/cF);
	glVertex3i(((LIST_Q[2][l]%(cH+1))-cH/2)/cF,(((LIST_Q[2][l]%((cH+1)*(cH+1)))/(cH+1))-cH/2)/cF,((LIST_Q[2][l]/((cH+1)*(cH+1)))-cH/2)/cF);
	glVertex3i(((LIST_Q[3][l]%(cH+1))-cH/2)/cF,(((LIST_Q[3][l]%((cH+1)*(cH+1)))/(cH+1))-cH/2)/cF,((LIST_Q[3][l]/((cH+1)*(cH+1)))-cH/2)/cF);
	glEnd();
	
}

glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_t()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_T[nump];
    
    glColor3f(1.,1.-pow(((tt-mintt)/(maxtt-mintt)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_phix()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_FX[nump];
    
    glColor3f(1.,1.-pow(((tt-minfx)/(maxfx-minfx)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_phiy()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_FY[nump];
    
    glColor3f(1.,1.-pow(((tt-minfy)/(maxfy-minfy)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_phiz()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_FZ[nump];
    
    glColor3f(1.,1.-pow(((tt-minfz)/(maxfz-minfz)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_ux()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_UX[nump];
    
    glColor3f(1.,1.-pow(((tt-minux)/(maxux-minux)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_uy()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_UY[nump];
    
    glColor3f(1.,1.-pow(((tt-minuy)/(maxuy-minuy)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_uz()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_UZ[nump];
    
    glColor3f(1.,1.-pow(((tt-minuz)/(maxuz-minuz)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_sxx()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_SXX[nump];
    
    glColor3f(1.,1.-pow(((tt-minsxx)/(maxsxx-minsxx)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_sxy()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_SXY[nump];
    
    glColor3f(1.,1.-pow(((tt-minsxy)/(maxsxy-minsxy)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_sxz()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_SXZ[nump];
    
    glColor3f(1.,1.-pow(((tt-minsxz)/(maxsxz-minsxz)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_syy()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_SYY[nump];
    
    glColor3f(1.,1.-pow(((tt-minsyy)/(maxsyy-minsyy)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_syz()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_SYZ[nump];
    
    glColor3f(1.,1.-pow(((tt-minsyz)/(maxsyz-minsyz)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);

    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

void afficher3_szz()
{
glClearColor(1.0,1.0,1.0,1.0);	
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(0,0,1.5*cH/cF,0,0,0,0,1,0);


//titre
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
if(dire==1){
const unsigned char myString [] = "Plan yz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==2){
const unsigned char myString [] = "Plan xz";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}else if(dire==3){
const unsigned char myString [] = "Plan xy";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 
}

for(int j=0;j<cH;j++){ 
 for(int i=0;i<cH;i++){

    int gamma,alpha,beta;
    
    if(dire==3){
    gamma = posi;  
    alpha = j;
    beta  = i;
    }
    else if(dire==2){
    gamma =j;
    alpha = posi;
    beta  = i;
    }
    else if(dire==1){
    gamma =j;
    alpha = i;
    beta  = posi;
    }    
   
    int nump  = gamma*cH*cH+alpha*cH+beta;
    R tt = LIST_SZZ[nump];
    
    glColor3f(1.,1.-pow(((tt-minszz)/(maxszz-minszz)),0.5),0.);  

    glBegin(GL_QUADS);

		glVertex2f((i-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j-cH/2.)/cF);
		glVertex2f((i+1-cH/2.)/cF,(j+1-cH/2.)/cF);
		glVertex2f((i-cH/2.)/cF,(j+1-cH/2.)/cF);
	
    glEnd(); 	//Pour les explications, lire le tutorial sur OGL et win
 
}  
}

// cadre map

	glColor3f(0.,0.,0.); 
	glBegin(GL_LINE_LOOP);
	glVertex2f((-cH/2.)/cF, (-cH/2.)/cF);     
	glVertex2f((cH/2.)/cF,  (-cH/2.)/cF);      
	glVertex2f((cH/2.)/cF,  (cH/2.)/cF);     
	glVertex2f((-cH/2.)/cF, (cH/2.)/cF);    
      	glEnd(); 

//info
glColor3f(0.,0.,0.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Position : " + to_string(posi+1) + "/" + to_string(cH) ;
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

      	
glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}



void refenetrerC(int w, int h)
{
glViewport(0, 0, (GLsizei) w, (GLsizei) h);
/*
modification des tailles
*/
/*du tampon d?affichage
*/
glMatrixMode(GL_PROJECTION);
/* pile courante = projection
*/
glLoadIdentity();
gluPerspective(
45,
float(w)/float(h),
0.1,
20*150
); 	//Pour les explications, lire le tutorial sur OGL et win
  	glMatrixMode(GL_MODELVIEW); 	//Optionnel


}



void refenetrer2(int w, int h)
{
glViewport(0, 0, (GLsizei) w, (GLsizei) h);
/*
modification des tailles
*/
/*du tampon d?affichage
*/
glMatrixMode(GL_PROJECTION);
/* pile courante = projection
*/
glLoadIdentity();
gluPerspective(
45,
float(w)/float(h),
0.1,
20*cH
); 	//Pour les explications, lire le tutorial sur OGL et win
  	glMatrixMode(GL_MODELVIEW); 	//Optionnel


}

void refenetrer3(int w, int h)
{
glViewport(0, 0, (GLsizei) w, (GLsizei) h);
/*
modification des tailles
*/
/*du tampon d?affichage
*/
glMatrixMode(GL_PROJECTION);
/* pile courante = projection
*/
glLoadIdentity();
gluPerspective(
45,
float(w)/float(h),
0.1,
40*cH
); 	//Pour les explications, lire le tutorial sur OGL et win
  	glMatrixMode(GL_MODELVIEW); 	//Optionnel


}

void clavier2(unsigned char touche, int x, int y)
{
	
	switch (touche)
	{
	case 27: {
	/* sortie si touche ESC */
              for(int ii=1;ii<=nb;ii++){       
	running[ii]= false;
	}
	glutLeaveMainLoop();
	break; 	}
    	case 43: {
	// touche +
	posi++;	
	if(posi>=cH) posi=cH-1;	
	
	break;
	}	
	case 45: {
	// touche -
	posi--;
	if(posi<0) posi=0;	
	
	break;
	}	
	case 'x': {
	dire=1;
	break;	
	}    
	case 'y': {
	dire=2;
	break;	
	}  
	case 'z': {
	dire=3;
	break;	
	}  

	}
}

void clavier3(unsigned char touche, int x, int y)
{

	switch (touche)
	{	
       
	case 27: {
	// sortie si touche ESC
        for(int ii=1;ii<=nb;ii++){       
	running[ii]= false;
	}
	glutLeaveMainLoop();
	break; }
	case 32: {
	// remise à 0 si space
	theta_eye=0;
	phi_eye=0;	
	rayon=2.5*((GLdouble) cH)/((GLdouble) cF);
	x_eye=0;
	y_eye=0;
	z_eye=rayon;
	x_g=280;
	y_g=280;
	x_d=280;
	y_d=280;
	break;
	}
	case 43: {
	// touche +
	posi++;
	
	if(posi>=cH) posi=cH-1;	
	
	break;
	}	
	case 45: {
	// touche -
	posi--;

	if(posi<0) posi=0;	
	
	break;
	}	

	case 81 : {
	// translation //z <0 si touche Q
        z_eye-=0.1*cH/cF;
	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	
	break;	  
	}
	case 113 : {
	// translation //z >0 si touche q
        z_eye+=0.1*cH/cF; 	
	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	break;	  
	}
	
	}
}

void processSpecialKeys(int key, int x, int y) {
  
        switch(key){
        case GLUT_KEY_LEFT: {
	// translation //X <0 si touche <-
        x_eye-=0.1*cH/cF;	

	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	
	break;	  
	}
        case GLUT_KEY_RIGHT : {
	// translation //X >0 si touche ->
        x_eye+=0.1*cH/cF;	

	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	
	break;	  
	}	
        case GLUT_KEY_DOWN : {
	// translation //Y <0 si touche v
        y_eye-=0.1*cH/cF;	

	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	
	break;	  
	}	
        case GLUT_KEY_UP : {
	// translation //Y >0 si touche ^
        y_eye+=0.1*cH/cF;	

	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	
	break;	  
	}	
	}
}


void gerer_souris2(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
cout<<"x_g : "<< x_g <<endl;
y_g     = 560-y;
cout<<"y_g : "<< y_g <<endl;  
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
cout<<"x_d : "<< x_d <<endl;
y_d     = 560-y;
cout<<"y_d : "<< y_d <<endl;    
}  
break;
}
}

void gerer_souris3(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
//cout<<"x_g : "<< x_g <<endl;
y_g     = 560-y;
//cout<<"y_g : "<< y_g <<endl;  
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
//cout<<"x_d : "<< x_d <<endl;
y_d     = 560-y;
//cout<<"y_d : "<< y_d <<endl;    
theta_eye+= (3.14159/4.*(double(y_d)-280.)/280.);
phi_eye+= (3.14159/4.*(double(x_d)-280.)/280.);  

x_eye = rayon*sin(3.14159/2-theta_eye)*sin(phi_eye);
y_eye = rayon*cos(3.14159/2-theta_eye);
z_eye = rayon*sin(3.14159/2-theta_eye)*cos(phi_eye);
  
}  
break;
}
}

void gerer_souris2_phix(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump  = alpha*cH+beta;
    R phitx = LIST_FX[nump];                
    cout<<" x :"<<(R(x_g-56)/448.)<<", y :"<<(R(y_g-56)/448.)<<", Flux x :"<<phitx<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris2_phiy(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump  = alpha*cH+beta;
    R phity = LIST_FY[nump];                
    cout<<" x :"<<(R(x_g-56)/448.)<<", y :"<<(R(y_g-56)/448.)<<", Flux y :"<<phity<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris2_t(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump  = alpha*cH+beta;
    R temp1 = LIST_T[nump];                
    cout<<" x :"<<(R(x_g-56)/448.)<<", y :"<<(R(y_g-56)/448.)<<", Temp :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris2_ux(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump  = alpha*cH+beta;
    R uu = LIST_UX[nump];                
    cout<<" x :"<<(R(x_g-56)/448.)<<", y :"<<(R(y_g-56)/448.)<<" Nx :"<<beta<<", Ny :"<<alpha<<", ux :"<<uu<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris2_uy(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump  = alpha*cH+beta;
    R uu = LIST_UY[nump];                
    cout<<" x :"<<(R(x_g-56)/448.)<<", y :"<<(R(y_g-56)/448.)<<", uy :"<<uu<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris2_sxx(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump  = alpha*cH+beta;
    R uu = LIST_SXX[nump];                
    cout<<" x :"<<(R(x_g-56)/448.)<<", y :"<<(R(y_g-56)/448.)<<", sxx :"<<uu<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris2_sxy(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump  = alpha*cH+beta;
    R uu = LIST_SXY[nump];                
    cout<<" x :"<<(R(x_g-56)/448.)<<", y :"<<(R(y_g-56)/448.)<<", sxy :"<<uu<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris2_syy(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump  = alpha*cH+beta;
    R uu = LIST_SYY[nump];                
    cout<<" x :"<<(R(x_g-56)/448.)<<", y :"<<(R(y_g-56)/448.)<<", syy :"<<uu<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_t(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_T[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", Temp :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_phix(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_FX[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", Phi :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_phiy(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_FY[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", Phi :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_phiz(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_FZ[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", Phi :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_sxx(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_SXX[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", ctr :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_sxy(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_SXY[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", ctr :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_sxz(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_SXZ[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", ctr :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_syy(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_SYY[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", ctr :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_syz(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_SYZ[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", ctr :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_szz(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_SZZ[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", ctr :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_ux(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_UX[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", dep :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_uy(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_UY[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", dep :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris3_uz(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
y_g     = 560-y;

    
    int beta = int((R(x_g-56)/448.)*cH);
    int alpha  = int((R(y_g-56)/448.)*cH);
    
    int nump;
    if(dire==3){
    nump = posi*cH*cH+alpha*cH+beta;
    }
    else if(dire==2){
    nump = alpha*cH*cH+posi*cH+beta;
    }   
    else if(dire==1){
    nump = alpha*cH*cH+beta*cH+posi;
    }       
    
    R temp1 = LIST_UZ[nump];                
    cout<<" x :"<<beta<<", y :"<<alpha<<", dep :"<<temp1<<endl;
    
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
y_d     = 560-y;
 
}  
break;
}
}

void gerer_souris_mouvement(int x, int y)
{
/*
position courante (x,y) de la souris
*/
}


void view_init2(int & argc, char ** argv){
glutInit(&argc, argv);
glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE);
glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);
 
InitGL2();
}

void view_init3(int & argc, char ** argv){
glutInit(&argc, argv);
glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE);
glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);
 
InitGL3();
}


void view_micro2_cont()
{
nb=nb+1;
running[nb] = true;
int nba=nb;
running[nba] = true;

glutCreateWindow("Microstructure 2D");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher2_cont);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris2);
glutMotionFunc(gerer_souris_mouvement);
while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();
}

void view_micro2_disk()
{
nb=nb+1;
running[nb] = true;
int nba=nb;
running[nba] = true;

glutCreateWindow("Microstructure 2D"); 

while(running[nba]){  
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher2_disk);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris2);
glutMotionFunc(gerer_souris_mouvement);
	
glutMainLoopEvent();
} 
glutDestroyWindow(nba+1);

//glutMainLoop();
}

void view_micro2_agr(){
nb=nb+1;
int nba=nb;
running[nba] = true;

glutCreateWindow("Microstructure 2D");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher2_agr);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris2);
glutMotionFunc(gerer_souris_mouvement);
while(running[nba]){  
	glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro2_t(){
nb=nb+1;

int nba=nb;
running[nba] = true;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_t);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de temperature");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher2_t);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris2_t);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro2_phix(){
nb=nb+1;

int nba=nb;
running[nba] = true;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_phix);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Flux direction x"); 
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher2_phix);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris2_phix);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro2_phiy(){
nb=nb+1;

int nba=nb;
running[nba] = true;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_phiy);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Flux direction y"); 
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher2_phiy);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris2_phiy);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro2_ux(){
nb=nb+1;

int nba=nb;
running[nba] = true;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_ux);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de deplacement ux"); 
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher2_ux);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris2_ux);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro2_uy(){
nb=nb+1;

int nba=nb;
running[nba] = true;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_uy);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de deplacement uy"); 
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher2_uy);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris2_uy);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro2_sxx(){
nb=nb+1;

int nba=nb;
running[nba] = true;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_sxx);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de contraintes sxx");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher2_sxx);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris2_sxx);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro2_sxy(){
nb=nb+1;

int nba=nb;
running[nba] = true;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_sxy);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de contraintes sxy");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher2_sxy);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris2_sxy);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro2_syy(){
nb=nb+1;

int nba=nb;
running[nba] = true;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_syy);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de contraintes syy");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher2_syy);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris2_syy);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}


void view_micro3_sphere()
{
nb=nb+1;
running[nb] = true;
int nba=nb;
running[nba] = true;

glutCreateWindow("Microstructure 3D"); 

while(running[nba]){  
glutReshapeFunc(refenetrer3);
glutDisplayFunc(afficher3_sphere);

glutKeyboardFunc(clavier3);
glutSpecialFunc(processSpecialKeys);
glutMouseFunc(gerer_souris3);
glutMotionFunc(gerer_souris_mouvement);
	
glutMainLoopEvent();
} 
glutDestroyWindow(nba+1);

//glutMainLoop();
}

void view_micro3_surf()
{
nb=nb+1;

int nba=nb;
running[nba] = true;

glutCreateWindow("Microstructure 3D"); 

while(running[nba]){  
glutReshapeFunc(refenetrer3);
glutDisplayFunc(afficher3_surf);

glutKeyboardFunc(clavier3);
glutSpecialFunc(processSpecialKeys);
glutMouseFunc(gerer_souris3);
glutMotionFunc(gerer_souris_mouvement);
	
glutMainLoopEvent();
} 
glutDestroyWindow(nba+1);

//glutMainLoop();
}


void view_micro3_t(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_t);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de temperature");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_t);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_t);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_phix(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_phix);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Flux direction x");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_phix);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_phix);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_phiy(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_phiy);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Flux direction y");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_phiy);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_phiy);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_phiz(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_phiz);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Flux direction z");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_phiz);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_phiz);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_sxx(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_sxx);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de contraintes sxx");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_sxx);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_sxx);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_sxy(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_sxy);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de contraintes sxy");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_sxy);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_sxy);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_sxz(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_sxz);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de contraintes sxz");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_sxz);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_sxz);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_syy(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_syy);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de contraintes syy");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_syy);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_syy);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_syz(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_syz);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de contraintes syz");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_syz);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_syz);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_szz(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_szz);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de contraintes szz");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_szz);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_szz);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_ux(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_ux);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de deplacement ux");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_ux);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_ux);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_uy(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_uy);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de deplacement uy");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_uy);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_uy);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}

void view_micro3_uz(int dire1, int posi1){
nb=nb+1;

int nba=nb;
running[nba] = true;

dire=dire1;
posi=posi1;

glutInitWindowSize(140,280);
glutInitWindowPosition(800,240);

glutCreateWindow("Colorbar");
glutReshapeFunc(refenetrerC);
glutDisplayFunc(afficherColorbar_uz);

nb=nb+1;
nba=nb;
running[nba] = true;

glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

glutCreateWindow("Champ de deplacement uz");
glutReshapeFunc(refenetrer2);
glutDisplayFunc(afficher3_uz);
  
glutKeyboardFunc(clavier2);
glutMouseFunc(gerer_souris3_uz);
glutMotionFunc(gerer_souris_mouvement);

while(running[nba]){  
glutMainLoopEvent();
} 
 glutDestroyWindow(nba+1);
//glutMainLoop();

}









