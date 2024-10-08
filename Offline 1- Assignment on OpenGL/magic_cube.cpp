#include <GL/glut.h>
#include<cmath>

GLint steps=16;
GLint rotateAngle=0;
struct point {
    GLfloat x, y, z;
};
struct point pos;   // position of the eye
struct point l;     // look/forward direction
struct point r;     // right direction
struct point u;     // up direction
GLdouble vertex1x=1,vertex1y=0,vertex1z=0,vertex2x=0,vertex2y=1,vertex2z=0,vertex3x=0,vertex3y=0,vertex3z=1;
GLdouble incrVal=(1-1.0/3)/steps;
GLdouble decrVal=(0-1.0/3)/steps;

GLdouble cylinderLen=sqrt((vertex1x-vertex2x)*(vertex1x-vertex2x)+(vertex1y-vertex2y)*(vertex1y-vertex2y)+(vertex1z-vertex2z)*(vertex1z-vertex2z));
//cylinderer circle xy plane e,so radius ber korte octahedronke clockwise(-45 degree) te ghrabo,x er basis e rotate so x fix thakbe,y,z change 
//hobe vertex1,2,3 sobgular jonno e
//vertex2 er basis e hishab kortesi tai just etake update korai enough
GLdouble cylinderCenterAngle = 70.5287794*M_PI/180;
GLdouble cylinderRadius=vertex2x/sin(cylinderCenterAngle/2);
//as vertex2 e asole rotate kortesina tai new instance nisi vertex2 er
GLdouble newVertex2y=(sqrt(2)/2) * vertex2y - (sqrt(2)/2) * vertex2z; //y' = cosθ * y - sinθ * z
//GLdouble newVertex2z= (-sqrt(2)/2) * vertex2y + (sqrt(2)/2) * vertex2z;//z' = sinθ * y + cosθ * z

GLdouble cylinderCentery=newVertex2y-vertex2x/tan(cylinderCenterAngle/2);//as 2d dhortesi,x=z=0, y er basis e (0,0,0) theke radius er distance y er basis e dhorbo
//sphere akar jnno 3 ta 2d array define korbo
const int rows = 33;
const int columns = 33;

double array2Dx[rows][columns];
double array2Dy[rows][columns];
double array2Dz[rows][columns];

GLdouble sphereRadius = sqrt(3)*vertex1y;
GLdouble sphereCenterx = vertex1x-vertex1y;


void drawTriangle()
{
    glBegin(GL_TRIANGLES);
    // glColor3f(1.0f, 0.0f, 0.0f); // Set triangle color to red
    glVertex3f(vertex1x,vertex1y,vertex1z); // Vertex 1
    glVertex3f(vertex2x,vertex2y,vertex2z); // Vertex 2
    glVertex3f(vertex3x,vertex3y,vertex3z); // Vertex 3
    glEnd();
}


void drawOctahedron()
{
    //glPushMatrix();
    glColor3f(0.0f,1.0f,1.0f);
     drawTriangle();
     //glPopMatrix();

    glPushMatrix();
    glRotatef(180.0f, 0.0f, 1.0f, 0.0f);
    
    drawTriangle();

    glColor3f(1.0f,0.0f,1.0f);
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    drawTriangle();
    glRotatef(180.0f, 0.0f, 1.0f, 0.0f);
    drawTriangle();
    glPopMatrix();

    glPushMatrix();
    glRotatef(180,1,0,0);
    glColor3f(0.0f,1.0f,1.0f);
     drawTriangle();
     //glPopMatrix();

    
    glRotatef(180.0f, 0.0f, 1.0f, 0.0f);
    
    drawTriangle();

    glColor3f(1.0f,0.0f,1.0f);
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    drawTriangle();
    glRotatef(180.0f, 0.0f, 1.0f, 0.0f);
    drawTriangle();
    glPopMatrix();

}

void drawCylinder(double height, double radius, int segments) {

    

    double tempx = radius, tempy = 0;
    double currx, curry;
    glBegin(GL_QUADS);
        for (int i = 1; i <= segments; i++) {
            double theta = i * 2.0 * M_PI / segments;
            currx = radius * cos(theta);
            curry = radius * sin(theta);

            GLfloat c = (2+cos(theta))/3;
            //glColor3f(c,c,c);
            glVertex3f(currx, curry, height/2);
            glVertex3f(currx, curry, -height/2);

            glVertex3f(tempx, tempy, -height/2);
            glVertex3f(tempx, tempy, height/2);

            tempx = currx;
            tempy = curry;
        }
    glEnd();
}

void drawFirstCylinder(){
     glPushMatrix();
    glRotatef(45,1,0,0);
    //glRotatef(45,-1,0,0);
    glTranslatef( 0,cylinderCentery,0);

    drawCylinder(cylinderLen,cylinderRadius,80);
    glPopMatrix();
}

void drawAllCylinder(){
    
    //glPushMatrix();
    glColor3f(1.0f,1.0f,0.0f);
    drawFirstCylinder();

    glPushMatrix();
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    
    drawFirstCylinder();
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    drawFirstCylinder();
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    drawFirstCylinder();
    glPopMatrix();

    glPushMatrix();
    glRotatef(180,1,0,0);
    drawFirstCylinder();
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    
    drawFirstCylinder();
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    drawFirstCylinder();
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    drawFirstCylinder();
    glPopMatrix();

    glPushMatrix();
    glRotatef(90,0,0,1);
    drawFirstCylinder();
    glRotatef(180.0f, 0.0f, 0.0f, 1.0f);
    drawFirstCylinder();
    glPopMatrix();
    
    glPushMatrix();
    glRotatef(180,0,1,0);
    glRotatef(90,0,0,1);
    drawFirstCylinder();
    glRotatef(180.0f, 0.0f, 0.0f, 1.0f);
    drawFirstCylinder();
    glPopMatrix();

}

void buildUnitPositiveX(int subdivision)
{
    const float DEG2RAD = acos(-1) / 180.0f;

    //std::vector<float> vertices;
    float n1[3];        // normal of longitudinal plane rotating along Y-axis
    float n2[3];        // normal of latitudinal plane rotating along Z-axis
    float v[3];         // direction vector intersecting 2 planes, n1 x n2
    float a1;           // longitudinal angle along Y-axis
    float a2;           // latitudinal angle along Z-axis

    // compute the number of vertices per row, 2^n + 1
    int pointsPerRow = (int)pow(2, subdivision) + 1;

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for(unsigned int i = 0; i < pointsPerRow; ++i)
    {
        // normal for latitudinal plane
        // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
        // therefore, it is rotating (0,1,0) vector by latitude angle a2
        a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
        n2[0] = -sin(a2);
        n2[1] = cos(a2);
        n2[2] = 0;

        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for(unsigned int j = 0; j < pointsPerRow; ++j)
        {
            // normal for longitudinal plane
            // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
            // therefore, it is rotating (0,0,-1) vector by longitude angle a1
            a1 = DEG2RAD * (-45.0f + 90.0f * j / (pointsPerRow - 1));
            n1[0] = -sin(a1);
            n1[1] = 0;
            n1[2] = -cos(a1);

            // find direction vector of intersected line, n1 x n2
            v[0] = n1[1] * n2[2] - n1[2] * n2[1];
            v[1] = n1[2] * n2[0] - n1[0] * n2[2];
            v[2] = n1[0] * n2[1] - n1[1] * n2[0];

            // normalize direction vector
            float scale = 1 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;

            
            array2Dx[i][j]=v[0];
            array2Dy[i][j]=v[1];
            array2Dz[i][j]=v[2];
        }
    }

   // return vertices;
}

void drawFirstSphere(){
    glPushMatrix();
    glTranslatef( sphereCenterx,0,0);
    glBegin(GL_QUADS);
    //glColor3f(1.0f, 0.0f, 0.0f);  // Set color to red
    for(int i=0;i<rows-1;i++){
        for(int j=0;j<columns-1;j++){
             glVertex3f(sphereRadius * array2Dx[i][j],sphereRadius* array2Dy[i][j],sphereRadius* array2Dz[i][j]);
             glVertex3f(sphereRadius *  array2Dx[i+1][j],sphereRadius *  array2Dy[i+1][j],sphereRadius *  array2Dz[i+1][j]);
             
             glVertex3f(sphereRadius *  array2Dx[i+1][j+1],sphereRadius *  array2Dy[i+1][j+1],sphereRadius *  array2Dz[i+1][j+1]);
             glVertex3f(sphereRadius *  array2Dx[i][j+1],sphereRadius *  array2Dy[i][j+1],sphereRadius *  array2Dz[i][j+1]);
             
        }
    }
    glEnd();
    glPopMatrix();

}

void drawAllSphere(){
    //xy plane er 4 ta sphere er 1/6th akbo
     glColor3f(0.0f, 1.0f, 0.0f);     // green
     drawFirstSphere();
     glPushMatrix();
     glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
     glColor3f(0.0f, 0.0f, 1.0f);     // Blue
     drawFirstSphere();
     glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
     glColor3f(0.0f, 1.0f, 0.0f);     // green
     drawFirstSphere();
     glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
     glColor3f(0.0f, 0.0f, 1.0f);     // Blue
     drawFirstSphere();
    glPopMatrix();
    //uporer 1/6th ar nicher 1/6th akbo
    glPushMatrix();
    glColor3f(1.0f, 0.0f, 0.0f);     // red
    glRotatef(90.0f, 0.0f, 0.0f, 1.0f);
    drawFirstSphere();
    glRotatef(180.0f, 0.0f, 0.0f, 1.0f);
    drawFirstSphere();
    glPopMatrix();

}

void keyboardListener(unsigned char key, int x, int y) {
    double rate = 0.05;
    double v = 0.25;
    //double lx = pos.x - l.x;
    //double lz = centerz - eyez;
    double s;
	switch(key){

		case '1':
			r.x = r.x*cos(rate)+l.x*sin(rate);
			r.y = r.y*cos(rate)+l.y*sin(rate);
			r.z = r.z*cos(rate)+l.z*sin(rate);

			l.x = l.x*cos(rate)-r.x*sin(rate);
			l.y = l.y*cos(rate)-r.y*sin(rate);
			l.z = l.z*cos(rate)-r.z*sin(rate);
			break;

        case '2':
			r.x = r.x*cos(-rate)+l.x*sin(-rate);
			r.y = r.y*cos(-rate)+l.y*sin(-rate);
			r.z = r.z*cos(-rate)+l.z*sin(-rate);

			l.x = l.x*cos(-rate)-r.x*sin(-rate);
			l.y = l.y*cos(-rate)-r.y*sin(-rate);
			l.z = l.z*cos(-rate)-r.z*sin(-rate);
			break;

        case '3':
			l.x = l.x*cos(rate)+u.x*sin(rate);
			l.y = l.y*cos(rate)+u.y*sin(rate);
			l.z = l.z*cos(rate)+u.z*sin(rate);

			u.x = u.x*cos(rate)-l.x*sin(rate);
			u.y = u.y*cos(rate)-l.y*sin(rate);
			u.z = u.z*cos(rate)-l.z*sin(rate);
			break;

        case '4':
			l.x = l.x*cos(-rate)+u.x*sin(-rate);
			l.y = l.y*cos(-rate)+u.y*sin(-rate);
			l.z = l.z*cos(-rate)+u.z*sin(-rate);

			u.x = u.x*cos(-rate)-l.x*sin(-rate);
			u.y = u.y*cos(-rate)-l.y*sin(-rate);
			u.z = u.z*cos(-rate)-l.z*sin(-rate);
			break;

        case '5':
			u.x = u.x*cos(rate)+r.x*sin(rate);
			u.y = u.y*cos(rate)+r.y*sin(rate);
			u.z = u.z*cos(rate)+r.z*sin(rate);

			r.x = r.x*cos(rate)-u.x*sin(rate);
			r.y = r.y*cos(rate)-u.y*sin(rate);
			r.z = r.z*cos(rate)-u.z*sin(rate);
			break;

        case '6':
			u.x = u.x*cos(-rate)+r.x*sin(-rate);
			u.y = u.y*cos(-rate)+r.y*sin(-rate);
			u.z = u.z*cos(-rate)+r.z*sin(-rate);

			r.x = r.x*cos(-rate)-u.x*sin(-rate);
			r.y = r.y*cos(-rate)-u.y*sin(-rate);
			r.z = r.z*cos(-rate)-u.z*sin(-rate);
			break;

        case 'a':
        
        rotateAngle-=10;
            rotateAngle%=360;
            
            break;
            

        case 'd':
            rotateAngle+=10;
            rotateAngle%=360;
            break;

    case ',':
        if(vertex1x - incrVal <1.0/3){
            vertex1x = vertex1y = vertex1z = vertex2x = vertex2y = vertex2z = vertex3x = vertex3y = vertex3z = 1.0/3; 
        }
        else{
            vertex1x=vertex1x-incrVal;
            vertex1y=vertex1y-decrVal;
            vertex1z=vertex1z-decrVal;
            vertex2y=vertex2y-incrVal;
            vertex2x=vertex2x-decrVal;
            vertex2z=vertex2z-decrVal;
            vertex3z=vertex3z-incrVal;
            vertex3x=vertex3x-decrVal;
            vertex3y=vertex3y-decrVal;
        }
        cylinderLen=sqrt((vertex1x-vertex2x)*(vertex1x-vertex2x)+(vertex1y-vertex2y)*(vertex1y-vertex2y)+(vertex1z-vertex2z)*(vertex1z-vertex2z));
        
        cylinderRadius=vertex2x/sin(cylinderCenterAngle/2);
        newVertex2y=(sqrt(2)/2) * vertex2y - (-sqrt(2)/2) * vertex2z;
        cylinderCentery=newVertex2y-vertex2x/tan(cylinderCenterAngle/2);
        sphereRadius = sqrt(3)*vertex1y;
        sphereCenterx = vertex1x-vertex1y;
        break;
    case '.':
        if(vertex1x + incrVal >1.0){
            vertex1x =  vertex2y =  vertex3z = 1.0; 
            vertex1y = vertex1z = vertex2x =vertex2z = vertex3x = vertex3y =0.0;
        }
        else{
            vertex1x=vertex1x+incrVal;
            vertex1y=vertex1y+decrVal;
            vertex1z=vertex1z+decrVal;
            vertex2y=vertex2y+incrVal;
            vertex2x=vertex2x+decrVal;
            vertex2z=vertex2z+decrVal;
            vertex3z=vertex3z+incrVal;
            vertex3x=vertex3x+decrVal;
            vertex3y=vertex3y+decrVal;
        }
        cylinderLen=sqrt((vertex1x-vertex2x)*(vertex1x-vertex2x)+(vertex1y-vertex2y)*(vertex1y-vertex2y)+(vertex1z-vertex2z)*(vertex1z-vertex2z));
        cylinderRadius=vertex2x/sin(cylinderCenterAngle/2);
        newVertex2y=(sqrt(2)/2) * vertex2y - (-sqrt(2)/2) * vertex2z;
        cylinderCentery=newVertex2y-vertex2x/tan(cylinderCenterAngle/2);
        sphereRadius = sqrt(3)*vertex1y;
        sphereCenterx = vertex1x-vertex1y;
        break;


    

    // Control exit
    case 27:    // ESC key
        exit(0);    // Exit window
        break;
    }
    glutPostRedisplay();    // Post a paint request to activate display()
}

void specialKeyListener(int key, int x,int y)
{
	switch(key){
		case GLUT_KEY_UP:		//down arrow key
			pos.x+=l.x;
			pos.y+=l.y;
			pos.z+=l.z;
			break;
		case GLUT_KEY_DOWN:		// up arrow key
			pos.x-=l.x;
			pos.y-=l.y;
			pos.z-=l.z;
			break;

		case GLUT_KEY_RIGHT:
			pos.x+=r.x;
			pos.y+=r.y;
			pos.z+=r.z;
			break;
		case GLUT_KEY_LEFT :
			pos.x-=r.x;
			pos.y-=r.y;
			pos.z-=r.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x+=u.x;
			pos.y+=u.y;
			pos.z+=u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
            pos.x-=u.x;
			pos.y-=u.y;
			pos.z-=u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
	glutPostRedisplay();
}



void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    //glTranslatef(0.0f, 0.0f, -5.0f); // Move the octahedron back along the z-axis

    // control viewing (or camera)
    gluLookAt(pos.x,pos.y,pos.z,
              pos.x+l.x,pos.y+l.y,pos.z+l.z,
              u.x,u.y,u.z);

    // drawFirstCylinder();
    glRotatef(rotateAngle, 0,1,0);
    drawAllCylinder();
    
    // glRotatef(-45,1,0,0);
drawOctahedron();
 drawAllSphere();
 

   


    glFlush();
}


/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0) height = 1;                // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    /*if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }*/
    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(45.0f, aspect, 0.1f, 50.0f);
}
// void reshape(int width, int height)
// {
//     glViewport(0, 0, width, height);
//     glMatrixMode(GL_PROJECTION);
//     glLoadIdentity();
//     gluPerspective(45.0, (double)width / (double)height, 1.0, 10.0);
//     glMatrixMode(GL_MODELVIEW);
// }

int main(int argc, char** argv)
{
    pos.x=2.5;pos.y=2.5;pos.z=2.5;

    l.x=-1.0/sqrt(3);l.y=-1.0/sqrt(3);l.z=-1.0/sqrt(3);
    u.x=-1.0/sqrt(6);u.y=2.0/sqrt(6);u.z=-1.0/sqrt(6);
    r.x=1.0/sqrt(2);r.y=0;r.z=-1.0/sqrt(2);
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Magic Cube");
    glEnable(GL_DEPTH_TEST);
    buildUnitPositiveX(5);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);

    glutMainLoop();

    return 0;
}
