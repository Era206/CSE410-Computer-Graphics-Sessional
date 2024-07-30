#include <GL/glut.h>
#include<cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "1805116_bitmap_image.hpp"
using namespace std;

GLint steps=16;
GLint rotateAngle=0;

struct Point{
    double x;
    double y;
    double z;

    Point() : x(0), y(0), z(0) {}
    Point(double xVal, double yVal, double zVal) : x(xVal), y(yVal), z(zVal) {}
};

Point operator+(const Point& p1, const Point& p2) {
    Point result(0.0,0.0,0.0);
    result.x = p1.x + p2.x;
    result.y = p1.y + p2.y;
    result.z = p1.z + p2.z;
    return result;
}

Point operator-(const Point& p1, const Point& p2) {
    Point result(0.0,0.0,0.0);
    result.x = p1.x - p2.x;
    result.y = p1.y - p2.y;
    result.z = p1.z - p2.z;
    return result;
}

Point operator*(const Point& p1, const Point& p2) {
    Point result(0.0,0.0,0.0);
    result.x = p1.x * p2.x;
    result.y = p1.y * p2.y;
    result.z = p1.z * p2.z;
    return result;
}

double dotProduct(const Point& p1, const Point& p2){
    double result=0.0;
    result=p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
    return result;
}

Point scalarProduct(const Point& p, double val){
    Point result(0.0,0.0,0.0);
    result.x=p.x*val;
    result.y=p.y*val;
    result.z=p.z*val;
    return result;
}

Point crossProduct(const Point& v1, const Point& v2) {
    Point result(0.0,0.0,0.0);
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}


Point normalize(const Point& p) {
    Point result(0.0,0.0,0.0);
    double length = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    result.x = p.x / length;
    result.y = p.y / length;
    result.z = p.z / length;
    return result;
}

struct Point pos;   // position of the eye
struct Point l;     // look/forward direction
struct Point r;     // right direction
struct Point u;     // up direction
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


//newly added parameters
double nearPlane, farPlane, fovY, aspectRatio;
int levelOfRecursion, numberOfPixels;
double checkerboardWidth;
double ambientCoefficient, diffuseCoefficient, reflectionCoefficient;
int numberOfObjects;
double windowHeight, windowWidth;
double numOflightSources;
double numOfSpotLightSources;

struct Color{
    double r;
    double g;
    double b;

    Color(double rVal, double gVal, double bVal) : r(rVal), g(gVal), b(bVal) {}
    Color() : r(0), g(0), b(0) {}
};

Color scalarProductColor(const Color& c, double val){
    Color result(0.0,0.0,0.0);
    result.r=c.r*val;
    result.g=c.g*val;
    result.b=c.b*val;
    return result;
}

Color operator+(const Color& c1, const Color& c2) {
    Color result(0.0,0.0,0.0);
    result.r = c1.r + c2.r;
    result.g = c1.g + c2.g;
    result.b = c1.b + c2.b;
    return result;
}


class IntersectedObject{
    public:
    Point intersectionPoint;
    Point normalVector;
    Color color;
    double t;
    double shininess;
    double ka;
    double kd;
    double ks;
    double kr;
   

    IntersectedObject() : intersectionPoint(0.0,0.0,0.0),normalVector(0.0,0.0,0.0),t(0.0),color(0.0,0.0,0.0), shininess(0.0), ka(0.0), kd(0.0),ks(0.0), kr(0.0)  {}
};




// sphere class
class Sphere {

public:
    Point center;
    double radius;
    Color color;
    
    double ka,kd,ks,kr;
    double shininess;
    Sphere(const Point& center, double radius, const Color& color, double ka, double kd, double ks, double kr, double shininess)
        : center(center), radius(radius), color(color), ka(ka), kd(kd), ks(ks), kr(kr), shininess(shininess){}

    void drawSphere() const{
        glPushMatrix();  // Push the current matrix
        // glColor3ub(color.r, color.g, color.b);  // Set color using RGB values
           glColor3ub(static_cast<unsigned char>(color.r * 255.0),
                   static_cast<unsigned char>(color.g * 255.0),
                   static_cast<unsigned char>(color.b * 255.0));

        glTranslated(center.x, center.y, center.z);  // Translate to sphere center
        glutSolidSphere(radius, 50, 50);  // Draw the solid sphere

        glPopMatrix();  // Pop the matrix
    }

    IntersectedObject intersect( Point position, Point dir) const {
        Point oc = operator-(position, center);
        double a = dotProduct(dir,dir);
        double b = 2.0 * dotProduct(oc,dir);
        double c = dotProduct(oc,oc) - radius * radius;
        double discriminant = b * b - 4 * a * c;
        double t;
        if (discriminant < 0) {
            t = -1.0;
        } else {
            double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
            double t2 = (-b + sqrt(discriminant)) / (2.0 * a);
            if(t1<t2){
                t=t1;
            }
            else{
                t=t2;
            }
            if(t<0){
                t=-1.0;
            }
        }

        IntersectedObject intersectedObject;
        intersectedObject.t=t;
        intersectedObject.intersectionPoint=operator+(position,scalarProduct(dir,t));
        intersectedObject.normalVector=normalize(operator-(intersectedObject.intersectionPoint,center));
        intersectedObject.color=color;
        return intersectedObject;


    }

};

vector<Sphere> spheres;

double determinant(const vector<vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    if (rows == 0 || cols == 0) {
        cerr << "Error: Matrix is empty." << endl;
        return 0.0;
    }

    double det = 0.0;

    if (rows == 1) {
        return matrix[0][0];
    }

    for (int i = 0; i < rows; ++i) {
        vector<vector<double>> subMatrix;
        for (int j = 1; j < rows; ++j) {
            vector<double> row;
            for (int k = 0; k < cols; ++k) {
                if (k != i) {
                    row.push_back(matrix[j][k]);
                }
            }
            subMatrix.push_back(row);
        }
        det += (i % 2 == 0 ? 1 : -1) * matrix[0][i] * determinant(subMatrix);
    }

    return det;
}

IntersectedObject intersecTionTriangle(Point position, Point dir, Point p1, Point p2, Point p3, Color color) {
    // double A[3][3] = {{p1.x - p2.x, p1.x - p3.x, dir.x},
    //                             {p1.y - p2.y, p1.y - p3.y, dir.y},
    //                             {p1.z - p2.z, p1.z - p3.z, dir.z}};
    IntersectedObject intersectedObject=IntersectedObject();
    double A[3][3], findBeta[3][3], findGamma[3][3], findT[3][3];

    // A declaratiom
    A[0][0] = p1.x - p2.x;
    A[0][1] = p1.x - p3.x;
    A[0][2] = dir.x;
    A[1][0] = p1.y - p2.y;
    A[1][1] = p1.y - p3.y;
    A[1][2] = dir.y;
    A[2][0] = p1.z - p2.z;
    A[2][1] = p1.z - p3.z;
    A[2][2] = dir.z;


    double detA = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                        A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                        A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
  
    if(detA==0){
        
        intersectedObject.t=-1.0;
        intersectedObject.intersectionPoint=Point(0.0,0.0,0.0);
        intersectedObject.normalVector=Point(0.0,0.0,0.0);
        intersectedObject.color=Color(0.0,0.0,0.0);
        return intersectedObject;
    }
    // vector<vector<double>> tempvect = {{5,3,7},{2,-5,8},{-6,4,9}};
   
    // vector<vector <double>> findBeta={{p1.x-position.x,p1.x-p3.x,dir.x},
    //                                 {p1.y-position.y,p1.y-p3.y,dir.y},
    //                                 {p1.z-position.z,p1.z-p3.z,dir.z}};
    // findBeta declaration
    findBeta[0][0] = p1.x - position.x;
    findBeta[0][1] = p1.x - p3.x;
    findBeta[0][2] = dir.x;
    findBeta[1][0] = p1.y - position.y;
    findBeta[1][1] = p1.y - p3.y;
    findBeta[1][2] = dir.y;
    findBeta[2][0] = p1.z - position.z;
    findBeta[2][1] = p1.z - p3.z;
    findBeta[2][2] = dir.z;

    double detBeta = findBeta[0][0] * (findBeta[1][1] * findBeta[2][2] - findBeta[1][2] * findBeta[2][1]) -
                        findBeta[0][1] * (findBeta[1][0] * findBeta[2][2] - findBeta[1][2] * findBeta[2][0]) +
                        findBeta[0][2] * (findBeta[1][0] * findBeta[2][1] - findBeta[1][1] * findBeta[2][0]);

    double beta=detBeta/detA;

    // vector<vector <double>> findGamma={{p1.x-p2.x,p1.x-position.x,dir.x},
    //                                 {p1.y-p2.y,p1.y-position.y,dir.y},
    //                                 {p1.z-p2.z,p1.z-position.z,dir.z}};

    // findGamma declaration
    findGamma[0][0] = p1.x - p2.x;
    findGamma[0][1] = p1.x - position.x;
    findGamma[0][2] = dir.x;
    findGamma[1][0] = p1.y - p2.y;
    findGamma[1][1] = p1.y - position.y;
    findGamma[1][2] = dir.y;
    findGamma[2][0] = p1.z - p2.z;
    findGamma[2][1] = p1.z - position.z;
    findGamma[2][2] = dir.z;

    double detGamma = findGamma[0][0] * (findGamma[1][1] * findGamma[2][2] - findGamma[1][2] * findGamma[2][1]) -
                        findGamma[0][1] * (findGamma[1][0] * findGamma[2][2] - findGamma[1][2] * findGamma[2][0]) +
                        findGamma[0][2] * (findGamma[1][0] * findGamma[2][1] - findGamma[1][1] * findGamma[2][0]);


    double gamma=detGamma/detA;

    // vector<vector <double>> findT={{p1.x-p2.x,p1.x-p3.x,p1.x-position.x},
    //                             {p1.y-p2.y,p1.y-p3.y,p1.y-position.y},
    //                             {p1.z-p2.z,p1.z-p3.z,p1.z-position.z}};

    // findT declaration
    findT[0][0] = p1.x - p2.x;
    findT[0][1] = p1.x - p3.x;
    findT[0][2] = p1.x - position.x;
    findT[1][0] = p1.y - p2.y;
    findT[1][1] = p1.y - p3.y;
    findT[1][2] = p1.y - position.y;
    findT[2][0] = p1.z - p2.z;
    findT[2][1] = p1.z - p3.z;
    findT[2][2] = p1.z - position.z;

    double detT = findT[0][0] * (findT[1][1] * findT[2][2] - findT[1][2] * findT[2][1]) -
                        findT[0][1] * (findT[1][0] * findT[2][2] - findT[1][2] * findT[2][0]) +
                        findT[0][2] * (findT[1][0] * findT[2][1] - findT[1][1] * findT[2][0]);

    double t=detT/detA;

    
    if(beta>=0 && gamma>=0 && beta+gamma<=1){
        
        intersectedObject.t=t;
        intersectedObject.intersectionPoint=operator+(position,scalarProduct(dir,t));
        intersectedObject.normalVector=normalize(crossProduct(operator-(p2,p1),operator-(p3,p1)));
        intersectedObject.color=color;
        
        
    }
    else{
       
        intersectedObject.t=-1.0;
        intersectedObject.intersectionPoint=Point(0.0,0.0,0.0);
        intersectedObject.normalVector=Point(0.0,0.0,0.0);
        intersectedObject.color=Color(0.0,0.0,0.0);
        // return intersectedObject;
    }
    return intersectedObject;
}




class Cube{
    public:
    Point lowerLeftPoint,lowerRightPoint,upperLeftPoint,upperRightPoint,lowerLeftPoint2,lowerRightPoint2,upperLeftPoint2,upperRightPoint2;
    Color color;
    double ka,kd,ks,kr;
    double shininess;
    Cube(const Point& lowerLeftPoint,const Point& lowerRightPoint,const Point& upperLeftPoint,const Point& upperRightPoint,const Point& lowerLeftPoint2,const Point& lowerRightPoint2,const Point& upperLeftPoint2,const Point& upperRightPoint2,const Color& color,double ka,double kd,double ks,double kr,double shininess)
        : lowerLeftPoint(lowerLeftPoint),lowerRightPoint(lowerRightPoint),upperLeftPoint(upperLeftPoint),upperRightPoint(upperRightPoint),lowerLeftPoint2(lowerLeftPoint2),lowerRightPoint2(lowerRightPoint2),upperLeftPoint2(upperLeftPoint2),upperRightPoint2(upperRightPoint2),color(color), ka(ka), kd(kd), ks(ks), kr(kr), shininess(shininess) {}
    void drawCube() const{
       
        glColor3ub(static_cast<unsigned char>(color.r * 255.0),
                   static_cast<unsigned char>(color.g * 255.0),
                   static_cast<unsigned char>(color.b * 255.0));

        
       


        glBegin(GL_QUADS);
        glVertex3f(lowerLeftPoint.x,lowerLeftPoint.y,lowerLeftPoint.z);
        glVertex3f(lowerRightPoint.x,lowerRightPoint.y,lowerRightPoint.z);
        glVertex3f(upperRightPoint.x,upperRightPoint.y,upperRightPoint.z);
        glVertex3f(upperLeftPoint.x,upperLeftPoint.y,upperLeftPoint.z);

        glVertex3f(lowerLeftPoint2.x,lowerLeftPoint2.y,lowerLeftPoint2.z);
        glVertex3f(lowerRightPoint2.x,lowerRightPoint2.y,lowerRightPoint2.z);
        glVertex3f(upperRightPoint2.x,upperRightPoint2.y,upperRightPoint2.z);
        glVertex3f(upperLeftPoint2.x,upperLeftPoint2.y,upperLeftPoint2.z);

        glVertex3f(lowerLeftPoint.x,lowerLeftPoint.y,lowerLeftPoint.z);
        glVertex3f(lowerRightPoint.x,lowerRightPoint.y,lowerRightPoint.z);
        glVertex3f(lowerRightPoint2.x,lowerRightPoint2.y,lowerRightPoint2.z);
        glVertex3f(lowerLeftPoint2.x,lowerLeftPoint2.y,lowerLeftPoint2.z);

        glVertex3f(upperLeftPoint.x,upperLeftPoint.y,upperLeftPoint.z);
        glVertex3f(upperRightPoint.x,upperRightPoint.y,upperRightPoint.z);
        glVertex3f(upperRightPoint2.x,upperRightPoint2.y,upperRightPoint2.z);
        glVertex3f(upperLeftPoint2.x,upperLeftPoint2.y,upperLeftPoint2.z);

        glVertex3f(lowerLeftPoint.x,lowerLeftPoint.y,lowerLeftPoint.z);
        glVertex3f(upperLeftPoint.x,upperLeftPoint.y,upperLeftPoint.z);
        glVertex3f(upperLeftPoint2.x,upperLeftPoint2.y,upperLeftPoint2.z);
        glVertex3f(lowerLeftPoint2.x,lowerLeftPoint2.y,lowerLeftPoint2.z);

        glVertex3f(lowerRightPoint.x,lowerRightPoint.y,lowerRightPoint.z);
        glVertex3f(upperRightPoint.x,upperRightPoint.y,upperRightPoint.z);
        glVertex3f(upperRightPoint2.x,upperRightPoint2.y,upperRightPoint2.z);
        glVertex3f(lowerRightPoint2.x,lowerRightPoint2.y,lowerRightPoint2.z);

        glEnd();
        // glPopMatrix();  // Pop the matrix
    }


    IntersectedObject intersect( Point position, Point dir) const {
        IntersectedObject intersectedObject=IntersectedObject();
        intersectedObject.t=farPlane;
        intersectedObject.intersectionPoint=Point(0.0,0.0,0.0);
        intersectedObject.normalVector=Point(0.0,0.0,0.0);
        intersectedObject.color=Color(0.0,0.0,0.0);

        
       
        // IntersectedObject intersectedObject1=intersecTionTriangle(position,dir,p1,p2,p3,color);
        //front face
        IntersectedObject intersectedObject1=intersecTionTriangle(position,dir,lowerLeftPoint,lowerRightPoint,upperLeftPoint,color);
    
        IntersectedObject intersectedObject2=intersecTionTriangle(position,dir,lowerRightPoint,upperRightPoint,upperLeftPoint,color);
        // up face
        IntersectedObject intersectedObject3=intersecTionTriangle(position,dir,upperLeftPoint,upperRightPoint,upperRightPoint2,color);
        IntersectedObject intersectedObject4=intersecTionTriangle(position,dir,upperLeftPoint2,upperLeftPoint,upperRightPoint2,color);
        // back face
        IntersectedObject intersectedObject5=intersecTionTriangle(position,dir,upperLeftPoint2,upperRightPoint2,lowerLeftPoint2,color);
        IntersectedObject intersectedObject6=intersecTionTriangle(position,dir,upperRightPoint2,lowerRightPoint2,lowerLeftPoint2,color);
        // down face
        IntersectedObject intersectedObject7=intersecTionTriangle(position,dir,lowerLeftPoint2,lowerRightPoint2,lowerRightPoint,color);
        IntersectedObject intersectedObject8=intersecTionTriangle(position,dir,lowerLeftPoint2,lowerRightPoint,lowerLeftPoint,color);
        // left face
        IntersectedObject intersectedObject9=intersecTionTriangle(position,dir,lowerLeftPoint,upperLeftPoint2,lowerLeftPoint2,color);
        IntersectedObject intersectedObject10=intersecTionTriangle(position,dir,lowerLeftPoint,upperLeftPoint,upperLeftPoint2,color);
        // right face
        IntersectedObject intersectedObject11=intersecTionTriangle(position,dir,lowerRightPoint2,upperRightPoint2,upperRightPoint,color);
        IntersectedObject intersectedObject12=intersecTionTriangle(position,dir,upperRightPoint,lowerRightPoint,lowerRightPoint2,color);
      

        if(intersectedObject1.t>0 && intersectedObject1.t<intersectedObject.t){
            intersectedObject=intersectedObject1;
       
        }
        if(intersectedObject2.t>0 && intersectedObject2.t<intersectedObject.t){
            intersectedObject=intersectedObject2;
        }
        
        if(intersectedObject3.t>0 && intersectedObject3.t<intersectedObject.t){
            intersectedObject=intersectedObject3;
        }
        if(intersectedObject4.t>0 && intersectedObject4.t<intersectedObject.t){
            intersectedObject=intersectedObject4;
        }
        if(intersectedObject5.t>0 && intersectedObject5.t<intersectedObject.t){
            intersectedObject=intersectedObject5;
        }
        if(intersectedObject6.t>0 && intersectedObject6.t<intersectedObject.t){
            intersectedObject=intersectedObject6;
        }
        if(intersectedObject7.t>0 && intersectedObject7.t<intersectedObject.t){
            intersectedObject=intersectedObject7;
            
        }
        if(intersectedObject8.t>0 && intersectedObject8.t<intersectedObject.t){
            intersectedObject=intersectedObject8;
           
        }
        if(intersectedObject9.t>0 && intersectedObject9.t<intersectedObject.t){
            intersectedObject=intersectedObject9;
            
        }
        if(intersectedObject10.t>0 && intersectedObject10.t<intersectedObject.t){
            intersectedObject=intersectedObject10;
            
        }
        if(intersectedObject11.t>0 && intersectedObject11.t<intersectedObject.t){
            intersectedObject=intersectedObject11;
           
        }
        if(intersectedObject12.t>0 && intersectedObject12.t<intersectedObject.t){
            intersectedObject=intersectedObject12;
           
        }
       
            
        return intersectedObject;
    }
};

vector <Cube> cubes;


class Pyramid{
    public:
    Point lowestPoint;
    double width,height;
    Color color;
    double ka,kd,ks,kr;
    double shininess;
    Point lowerRightPoint,upperLeftPoint,upperRightPoint,lowerLeftPoint,topPoint;
    Pyramid(const Point& lowestPoint,double width,double height,const Color& color, double ka, double kd, double ks, double kr, double shininess)
        : lowestPoint(lowestPoint),width(width),height(height),color(color),ka(ka), kd(kd), ks(ks), kr(kr), shininess(shininess) {}

    // need to implement a function to draw pyrimid, so I have e lowest point, that is the middle point of the base of the pyrimid, I need to find the other 4 points of the base and the top point of the pyrimid
    // I can find the other 4 points of the base by using the width and height of the pyrimid, I can find the top point by using the lowest point and the height of the pyrimid
    // then I can draw the pyrimid using the 5 points, now let's implement the function
    void calculatePoints(){
         lowerRightPoint.x=lowestPoint.x+width/2;
        lowerRightPoint.y=lowestPoint.y+width/2;
        lowerRightPoint.z=lowestPoint.z;

        lowerLeftPoint.x=lowestPoint.x+width/2;
        lowerLeftPoint.y=lowestPoint.y-width/2;
        lowerLeftPoint.z=lowestPoint.z;

        upperRightPoint.x=lowestPoint.x-width/2;
        upperRightPoint.y=lowestPoint.y+width/2;
        upperRightPoint.z=lowestPoint.z;

        upperLeftPoint.x=lowestPoint.x-width/2;
        upperLeftPoint.y=lowestPoint.y-width/2;
        upperLeftPoint.z=lowestPoint.z;

        topPoint.x=lowestPoint.x;
        topPoint.y=lowestPoint.y;
        topPoint.z=lowestPoint.z+height;
    }
    void drawPyramid() const{
        
        // glPushMatrix();  // Push the current matrix
         glColor3ub(static_cast<unsigned char>(color.r * 255.0),
                   static_cast<unsigned char>(color.g * 255.0),
                   static_cast<unsigned char>(color.b * 255.0));
        glBegin(GL_QUADS);
        glVertex3f(lowerLeftPoint.x,lowerLeftPoint.y,lowerLeftPoint.z);
        glVertex3f(lowerRightPoint.x,lowerRightPoint.y,lowerRightPoint.z);
        glVertex3f(upperRightPoint.x,upperRightPoint.y,upperRightPoint.z);
        glVertex3f(upperLeftPoint.x,upperLeftPoint.y,upperLeftPoint.z);

        glEnd();
        glBegin(GL_TRIANGLES);
        glVertex3f(lowerLeftPoint.x,lowerLeftPoint.y,lowerLeftPoint.z);
        glVertex3f(lowerRightPoint.x,lowerRightPoint.y,lowerRightPoint.z);
        glVertex3f(topPoint.x,topPoint.y,topPoint.z);

        glVertex3f(lowerRightPoint.x,lowerRightPoint.y,lowerRightPoint.z);
        glVertex3f(upperRightPoint.x,upperRightPoint.y,upperRightPoint.z);
        glVertex3f(topPoint.x,topPoint.y,topPoint.z);

        glVertex3f(upperRightPoint.x,upperRightPoint.y,upperRightPoint.z);
        glVertex3f(upperLeftPoint.x,upperLeftPoint.y,upperLeftPoint.z);
        glVertex3f(topPoint.x,topPoint.y,topPoint.z);

        glVertex3f(upperLeftPoint.x,upperLeftPoint.y,upperLeftPoint.z);
        glVertex3f(lowerLeftPoint.x,lowerLeftPoint.y,lowerLeftPoint.z);
        glVertex3f(topPoint.x,topPoint.y,topPoint.z);


        glEnd();

    }

    IntersectedObject intersect(Point position, Point dir) const{
        Point p1,p2,p3,p4,p5;

        p1=lowerLeftPoint;
        p2=lowerRightPoint;
        p3=upperRightPoint;
        p4=upperLeftPoint;
        p5=topPoint;
    

        IntersectedObject intersectedObject;
        intersectedObject.t=farPlane;
        intersectedObject.intersectionPoint=Point(0.0,0.0,0.0);
        intersectedObject.normalVector=Point(0.0,0.0,0.0);
        intersectedObject.color=Color(0.0,0.0,0.0);

        IntersectedObject intersectedObject1=intersecTionTriangle(position,dir,p1,p2,p5,color);
        IntersectedObject intersectedObject2=intersecTionTriangle(position,dir,p2,p3,p5,color);
        IntersectedObject intersectedObject3=intersecTionTriangle(position,dir,p3,p4,p5,color);
        IntersectedObject intersectedObject4=intersecTionTriangle(position,dir,p4,p1,p5,color);
        IntersectedObject intersectedObject5=intersecTionTriangle(position,dir,p1,p3,p2,color);
        IntersectedObject intersectedObject6=intersecTionTriangle(position,dir,p1,p4,p3,color);

        if(intersectedObject1.t>0){
            intersectedObject=intersectedObject1;
        }
        if(intersectedObject2.t>0 && intersectedObject2.t<intersectedObject.t){
            intersectedObject=intersectedObject2;
        }
        if(intersectedObject3.t>0 && intersectedObject3.t<intersectedObject.t){
            intersectedObject=intersectedObject3;
        }
        if(intersectedObject4.t>0 && intersectedObject4.t<intersectedObject.t){
            intersectedObject=intersectedObject4;
        }
        if(intersectedObject5.t>0 && intersectedObject5.t<intersectedObject.t){
            intersectedObject=intersectedObject5;
        }
        if(intersectedObject6.t>0 && intersectedObject6.t<intersectedObject.t){
            intersectedObject=intersectedObject6;
        }
        return intersectedObject;
    }

    

};

vector <Pyramid> pyramids;

class pointLight{
    public:
    Point position;
    double falloff;
    
    pointLight(const Point& position, double falloff)
        : position(position), falloff(falloff) {}
        
    pointLight() {}

    void drawPointLight() const{
        glPushMatrix();
        glTranslatef(position.x,position.y,position.z);
        glColor3f(1.0f, 1.0f, 1.0f);
        glutSolidSphere(30, 60, 60);
        glPopMatrix();
    }
};

vector <pointLight> pointLights;

class spotLight{
    public:
    Point position;
    Point direction;
    double cutoffAngle;
    double falloff;
    
    spotLight(const Point& position, const Point& direction, double cutoffAngle, double falloff)
        : position(position), direction(direction), cutoffAngle(cutoffAngle), falloff(falloff) {}

    spotLight() {}
    
    void drawSpotLight() const{
        glPushMatrix();
        glTranslatef(position.x,position.y,position.z);
        glColor3f(1.0f, 1.0f, 1.0f);
        glutSolidSphere(30, 60, 60);
        glPopMatrix();
    }
};

vector <spotLight> spotLights;


// Helper function to parse a line and extract space-separated values
vector<double> parseLine(const  string& line) {
     vector<double> values;
     istringstream iss(line);
    double value;
    while (iss >> value) {
        values.push_back(value);
    }
    return values;
}

void readFromFile(const  string& filename){
    
    ifstream infile(filename);

    // Read fixed parameters
    for (int i = 0; i < 6; ++i) {
         string line;
         getline(infile, line);
        if(!line.empty()){
            vector<double> values = parseLine(line);
            if(i==0){
                nearPlane=values[0];
                farPlane=values[1];
                fovY=values[2];
                aspectRatio=values[3];
            }
            if(i==1){
                levelOfRecursion=values[0];
            }
            if(i==2){
                numberOfPixels=values[0];
            }
            if(i==3){
                checkerboardWidth=values[0];
            }
            if(i==4){
                ambientCoefficient=values[0];
                diffuseCoefficient=values[1];
                // specularCoefficient=values[2];
                reflectionCoefficient=values[2];
            }
            if(i==5){
                
                numberOfObjects=values[0];
                
                // while (getline(infile, line)) {
                    // vector<double> values = parseLine(line);
                    // Do something with the values
                getline(infile, line);
                for(int n=0;n<numberOfObjects;n++){
                   
                    string line1;
                    getline(infile, line1);
                    while(line1.empty()){
                        getline(infile, line1);
                    }
                    if(line1=="cube"){
                        
                        Point cubeLowerLeftPoint(0.0,0.0,0.0);
                        double cubeSide=0.0;
                        Color cubeColor(0.0,0.0,0.0);
                        double cubeAmbientCoefficient,cubeDiffuseCoefficient,cubeSpecularCoefficient, cubeReflectionCoefficient;
                        double cubeShininess;
                        for(int j=0;j<5;j++){
                            string lines;
                            getline(infile, lines);
                            vector<double> values = parseLine(lines);
                            if(j==0){
                                cubeLowerLeftPoint.x=values[0];
                                cubeLowerLeftPoint.y=values[1];
                                cubeLowerLeftPoint.z=values[2];
                            }
                            if(j==1){
                                cubeSide=values[0];
                            }
                            if(j==2){
                                cubeColor.r=values[0];
                                cubeColor.g=values[1];
                                cubeColor.b=values[2];
                            }
                            if(j==3){
                                cubeAmbientCoefficient=values[0];
                                cubeDiffuseCoefficient=values[1];
                                cubeSpecularCoefficient=values[2];
                                cubeReflectionCoefficient=values[3];
                            }
                            if(j==4){
                                cubeShininess=values[0];
                            }

                        }
                        Point cubeLowerrightPoint(0.0,0.0,0.0);
                        cubeLowerrightPoint.x=cubeLowerLeftPoint.x;
                        cubeLowerrightPoint.y=cubeLowerLeftPoint.y+cubeSide;
                        cubeLowerrightPoint.z=cubeLowerLeftPoint.z;
                        Point cubeUpperleftPoint(0.0,0.0,0.0);
                        cubeUpperleftPoint.x=cubeLowerLeftPoint.x;
                        cubeUpperleftPoint.y=cubeLowerLeftPoint.y;
                        cubeUpperleftPoint.z=cubeLowerLeftPoint.z+cubeSide;
                        Point cubeUpperRightPoint(0.0,0.0,0.0);
                        cubeUpperRightPoint.x=cubeLowerLeftPoint.x;
                        cubeUpperRightPoint.y=cubeLowerLeftPoint.y+cubeSide;
                        cubeUpperRightPoint.z=cubeLowerLeftPoint.z+cubeSide;
                        Point cubeLowerLeftPoint2(0.0,0.0,0.0);
                        cubeLowerLeftPoint2.x=cubeLowerLeftPoint.x-cubeSide;
                        cubeLowerLeftPoint2.y=cubeLowerLeftPoint.y;
                        cubeLowerLeftPoint2.z=cubeLowerLeftPoint.z;
                        Point cubeLowerrightPoint2(0.0,0.0,0.0);
                        cubeLowerrightPoint2.x=cubeLowerLeftPoint.x-cubeSide;
                        cubeLowerrightPoint2.y=cubeLowerLeftPoint.y+cubeSide;
                        cubeLowerrightPoint2.z=cubeLowerLeftPoint.z;
                        Point cubeUpperleftPoint2(0.0,0.0,0.0);
                        cubeUpperleftPoint2.x=cubeLowerLeftPoint.x-cubeSide;
                        cubeUpperleftPoint2.y=cubeLowerLeftPoint.y;
                        cubeUpperleftPoint2.z=cubeLowerLeftPoint.z+cubeSide;
                        Point cubeUpperRightPoint2(0.0,0.0,0.0);
                        cubeUpperRightPoint2.x=cubeLowerLeftPoint.x-cubeSide;
                        cubeUpperRightPoint2.y=cubeLowerLeftPoint.y+cubeSide;
                        cubeUpperRightPoint2.z=cubeLowerLeftPoint.z+cubeSide;

                        Cube cube(cubeLowerLeftPoint,cubeLowerrightPoint,cubeUpperleftPoint,cubeUpperRightPoint,cubeLowerLeftPoint2,cubeLowerrightPoint2,cubeUpperleftPoint2,cubeUpperRightPoint2,cubeColor, cubeAmbientCoefficient, cubeDiffuseCoefficient, cubeSpecularCoefficient, cubeReflectionCoefficient, cubeShininess);
                        cubes.push_back(cube);
                    }
                    else if(line1=="sphere"){
                       
                        Point sphereCenter(0.0,0.0,0.0);
                        double sphereRadius;
                        Color sphereColor(0.0,0.0,0.0);
                        double sphereAmbientCoefficient,sphereDiffuseCoefficient,sphereSpecularCoefficient,sphereReflectionCoefficient;
                        double sphereShininess;
                        for(int k=0;k<6;k++){
                            
                            string lines;
                            getline(infile, lines);
                            vector<double> values = parseLine(lines);
                            if(k==0){
                                sphereCenter.x=values[0];
                                sphereCenter.y=values[1];
                                sphereCenter.z=values[2];
                            }
                            if(k==1){
                                sphereRadius=values[0];
                            }
                            if(k==2){
                                sphereColor.r=values[0];
                                sphereColor.g=values[1];
                                sphereColor.b=values[2];
                            }
                            if(k==3){
                                sphereAmbientCoefficient=values[0];
                                sphereDiffuseCoefficient=values[1];
                                sphereSpecularCoefficient=values[2];
                                sphereReflectionCoefficient=values[3];
                            }
                            if(k==4){
                                sphereShininess=values[0];
                            }
                        }
                        Sphere sphere(sphereCenter, sphereRadius, sphereColor, sphereAmbientCoefficient, sphereDiffuseCoefficient, sphereSpecularCoefficient, sphereReflectionCoefficient, sphereShininess );
                        spheres.push_back(sphere);

                    }
                    else{
                       
                        Point pyramidLowestPoint(0.0,0.0,0.0);
                        double pyramidWidth,pyramidHeight;
                        Color pyramidColor(0.0,0.0,0.0);
                        double pyramidAmbientCoefficient,pyramidDiffuseCoefficient,pyramidSpecularCoefficient,pyramidReflectionCoefficient;
                        double pyramidShininess;
                        for(int l=0;l<6;l++){
                            
                            string lines;
                            getline(infile, lines);
                            vector<double> values = parseLine(lines);
                            if(l==0){
                                pyramidLowestPoint.x=values[0];
                                pyramidLowestPoint.y=values[1];
                                pyramidLowestPoint.z=values[2];
                                
                            }
                            if(l==1){
                                pyramidWidth=values[0];
                                pyramidHeight=values[1];
                            }
                            if(l==2){
                                pyramidColor.r=values[0];
                                pyramidColor.g=values[1];
                                pyramidColor.b=values[2];
                            }
                            if(l==3){
                                pyramidAmbientCoefficient=values[0];
                                pyramidDiffuseCoefficient=values[1];
                                pyramidSpecularCoefficient=values[2];
                                pyramidReflectionCoefficient=values[3];
                            }
                            if(l==4){
                                pyramidShininess=values[0];
                            }
                        }
                        pyramidLowestPoint.x=pyramidLowestPoint.x-pyramidWidth/2;
                        pyramidLowestPoint.y=pyramidLowestPoint.y+pyramidWidth/2;
                    


                        Pyramid pyramid(pyramidLowestPoint,pyramidWidth,pyramidHeight,pyramidColor, pyramidAmbientCoefficient, pyramidDiffuseCoefficient, pyramidSpecularCoefficient, pyramidReflectionCoefficient, pyramidShininess);
                        pyramid.calculatePoints();
                        pyramids.push_back(pyramid);
                    }
            }


            }
            
        }
        else{
            i=i-1;
        }  
        
    }
    string line;
    getline(infile, line);
    while(line.empty()){
        getline(infile, line);
    }
    //    if(!line.empty()){
         vector<double> values = parseLine(line);
        numOflightSources=values[0];
        for(int j=0;j<numOflightSources;j++){
            getline(infile, line);
            vector<double> values = parseLine(line);
            Point point(values[0],values[1],values[2]);
            getline(infile, line);
            values = parseLine(line);
            pointLight light(point, values[0]);
            pointLights.push_back(light);
            getline(infile, line);
        }
    //    }
    getline(infile, line);
    while(line.empty()){
        getline(infile, line);
    }
        values = parseLine(line);
       
        numOfSpotLightSources=values[0];
        for(int j=0;j<numOfSpotLightSources;j++){
            getline(infile, line);
            vector<double> values = parseLine(line);
            Point point(values[0],values[1],values[2]);

            getline(infile, line);
            values = parseLine(line);
            double fallOffAngle=values[0];

            getline(infile, line);
            values = parseLine(line);
            Point direction(values[0],values[1],values[2]);
            
            // getline(infile, line);
            // values = parseLine(line);

            spotLight light(point,direction, values[3],fallOffAngle);
            spotLights.push_back(light);
            getline(infile, line);
        }

     while (getline(infile, line)) {
        if (line == "input explanation") {
            break; 
        }
    infile.close();
     }


}





void drawXZPlane(float width) {
    glBegin(GL_QUADS);
    int x = 0;
    for (int i = -500; i < 500; i++) {
        for (int k = -500; k < 500; k++) {
            if ((i + k) % 2 == 0) {
                x = 1;
                glColor3f(1.0f, 1.0f, 1.0f);  // Set color to white
            } else {
                x = 0;
                glColor3f(0.0f, 0.0f, 0.0f);  // Set color to black
            }
            glVertex3f(i * width, 0.0f, k * width);
            glVertex3f(i * width, 0.0f, (k + 1) * width);
            glVertex3f((i + 1) * width, 0.0f, (k + 1) * width);
            glVertex3f((i + 1) * width, 0.0f, k * width);
        }
    }
    glEnd();
}

void drawXYPlane(double width){
    glBegin(GL_QUADS);
    int x=0;
    for(int i=-500;i<500;i=i+1){
        for(int j=-500;j<500;j=j+1){
            if((i+j)%2==0){
                x=1;
                glColor3f(1.0f, 1.0f, 1.0f);  // Set color to white
            }
            else{
                
                glColor3f(0.0f, 0.0f, 0.0f);  // Set color to black
            }
           
            // draw checkerboard in XY plane
            glVertex3f(i*width,j*width,0);
            glVertex3f(i*width,(j+1.0)*width,0);
            glVertex3f((i+1.0)*width,(j+1.0)*width,0);
            glVertex3f((i+1.0)*width,j*width,0);
        }
    }           
    glEnd();
}




void drawAxis()
{
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);  // Set color to red
    glVertex3f(0.0f,0.0f,0.0f);
    glVertex3f(100.0f,0.0f,0.0f);
    glColor3f(0.0f, 1.0f, 0.0f);  // Set color to green
    glVertex3f(0.0f,0.0f,0.0f);
    glVertex3f(0.0f,100.0f,0.0f);
    glColor3f(0.0f, 0.0f, 1.0f);  // Set color to blue
    glVertex3f(0.0f,0.0f,0.0f);
    glVertex3f(0.0f,0.0f,100.0f);
    glEnd();
}

void drawTriangle()
{
    glBegin(GL_TRIANGLES);
    // glColor3f(1.0f, 0.0f, 0.0f); // Set triangle color to red
    glVertex3f(vertex1x,vertex1y,vertex1z); // Vertex 1
    glVertex3f(vertex2x,vertex2y,vertex2z); // Vertex 2
    glVertex3f(vertex3x,vertex3y,vertex3z); // Vertex 3
    glEnd();
}

IntersectedObject intersectCheckerBoard(Point position, Point direction){
    IntersectedObject answer=IntersectedObject();
    answer.t=-1.0;
    if(direction.z==0){
        return answer;
    }
    double t=-position.z/direction.z;
    if(t<0){
        return answer;
    }
    Point intersectionPoint=operator+(position,scalarProduct(direction, t));
    int baseX=(int)(floor(intersectionPoint.x/checkerboardWidth));
    int baseY=(int)(floor(intersectionPoint.y/checkerboardWidth));
    if((baseX+baseY)%2==0){
        answer.color=Color(1.0,1.0,1.0);
    }
    else{
        answer.color=Color(0.0,0.0,0.0);
    }
    answer.t=t;
    answer.intersectionPoint=intersectionPoint;
    if(position.z<0){
        answer.normalVector=Point(0.0,0.0,-1.0);
    }
    else{
        answer.normalVector=Point(0.0,0.0,1.0);
    }
    return answer;
}

double pointBuffer[1000][1000][3];

IntersectedObject doesIlluminate(Point pixelPoint, Point dirVector){
    IntersectedObject intersectedObject;
    intersectedObject.t=farPlane;
 
    Color color(0.0,0.0,0.0);
    for(int k=0;k<spheres.size();k++){
        IntersectedObject intersectedObject1=spheres[k].intersect(pixelPoint,dirVector);
        if(intersectedObject1.t>0 && intersectedObject1.t<intersectedObject.t){
            // nearest=k;
            intersectedObject=intersectedObject1;
            color=spheres[k].color;
            intersectedObject.color = color;
        }
    }
    // printf("Here 2\n");
    for(int k=0;k<cubes.size();k++){
        IntersectedObject intersectedObject1=cubes[k].intersect(pixelPoint,dirVector);
        if(intersectedObject1.t>0 && intersectedObject1.t<intersectedObject.t){
            // nearest=k;
            intersectedObject=intersectedObject1;
            color=cubes[k].color;
            intersectedObject.color = color;
           
        }
    }
    // printf("Here 3\n");
    for(int k=0;k<pyramids.size();k++){
        IntersectedObject intersectedObject1=pyramids[k].intersect(pixelPoint,dirVector);
        if(intersectedObject1.t>0 && intersectedObject1.t<intersectedObject.t){
            // nearest=k;
            intersectedObject=intersectedObject1;
            color=pyramids[k].color;
            intersectedObject.color = color;
        }
    }


    IntersectedObject intersectedObject2=intersectCheckerBoard(pixelPoint,dirVector);
    if(intersectedObject2.t>0 && intersectedObject2.t<intersectedObject.t){
        // nearest=-2;
        intersectedObject=intersectedObject2;
        color=intersectedObject2.color;
        intersectedObject.color = color;
    }

    return intersectedObject;
}



IntersectedObject rayTrace(Point pixelPoint, Point dirVector, int levelOfRecursion){
    IntersectedObject intersectedObject;
    intersectedObject.t=farPlane;
    // printf("Here\n");
    // tMin=farPlane;
    // nearest=-1;
    Color color(0.0,0.0,0.0);
    for(int k=0;k<spheres.size();k++){
        IntersectedObject intersectedObject1=spheres[k].intersect(pixelPoint,dirVector);
        if(intersectedObject1.t>0 && intersectedObject1.t<intersectedObject.t){
            // nearest=k;
            intersectedObject=intersectedObject1;
            color=spheres[k].color;
            intersectedObject.color = color;
            intersectedObject.shininess=spheres[k].shininess;
            intersectedObject.ka=spheres[k].ka;
            intersectedObject.kd=spheres[k].kd;
            intersectedObject.ks=spheres[k].ks;
            intersectedObject.kr=spheres[k].kr;
            
        }
    }
    // printf("Here 2\n");
    for(int k=0;k<cubes.size();k++){
        IntersectedObject intersectedObject1=cubes[k].intersect(pixelPoint,dirVector);
        if(intersectedObject1.t>0 && intersectedObject1.t<intersectedObject.t){
            // nearest=k;
            intersectedObject=intersectedObject1;
            color=cubes[k].color;
            intersectedObject.color = color;
            intersectedObject.shininess=cubes[k].shininess;
            intersectedObject.ka=cubes[k].ka;
            intersectedObject.kd=cubes[k].kd;
            intersectedObject.ks=cubes[k].ks;
            intersectedObject.kr=cubes[k].kr;
            
        }
    }
    // printf("Here 3\n");
    for(int k=0;k<pyramids.size();k++){
        IntersectedObject intersectedObject1=pyramids[k].intersect(pixelPoint,dirVector);
        if(intersectedObject1.t>0 && intersectedObject1.t<intersectedObject.t){
            // nearest=k;
            intersectedObject=intersectedObject1;
            color=pyramids[k].color;
            intersectedObject.color = color;
            intersectedObject.shininess=pyramids[k].shininess;
            intersectedObject.ka=pyramids[k].ka;
            intersectedObject.kd=pyramids[k].kd;
            intersectedObject.ks=pyramids[k].ks;
            intersectedObject.kr=pyramids[k].kr;

        }
    }
    // printf("Here 4\n");
    // printf("color: %lf %lf %lf\n",color.r,color.g,color.b);

    IntersectedObject intersectedObject2=intersectCheckerBoard(pixelPoint,dirVector);
    if(intersectedObject2.t>0 && intersectedObject2.t<intersectedObject.t){
        // nearest=-2;
        intersectedObject=intersectedObject2;
        color=intersectedObject2.color;
        intersectedObject.color = color;
        intersectedObject.shininess=0.0;
        intersectedObject.ka=ambientCoefficient;
        intersectedObject.kd=diffuseCoefficient;
        intersectedObject.ks=0.0;
        intersectedObject.kr=reflectionCoefficient;
    }

    if(intersectedObject.t == farPlane){
        intersectedObject.color=Color(0.0,0.0,0.0);
        intersectedObject.t = -1;
        return intersectedObject;
    }

    double lambert=0.0;
    double phong=0.0;
    for(int i=0;i<pointLights.size();i++){
        Point pointToLightVector=operator-(pointLights[i].position,intersectedObject.intersectionPoint);
        pointToLightVector=normalize(pointToLightVector);
        
        Point lightToPointVector=operator-(intersectedObject.intersectionPoint,pointLights[i].position);
        lightToPointVector=normalize(lightToPointVector);

        Point normalVector=intersectedObject.normalVector;
        normalVector=normalize(normalVector);

        Point lightReflectedVector=operator-(lightToPointVector,scalarProduct(normalVector,2.0*dotProduct(normalVector,lightToPointVector)));

        double t1=(pointLights[i].position.x-intersectedObject.intersectionPoint.x)/pointToLightVector.x;

        IntersectedObject intersectedObject3=doesIlluminate(pointLights[i].position,lightToPointVector);

        double t2=intersectedObject3.t;

        if(((t2-t1)*(t2-t1))>0.0000001){
            continue;
        }

        double scalingFactor = exp(-t1*t1*pointLights[i].falloff);
        lambert+=scalingFactor*dotProduct(normalVector,pointToLightVector);
        
        Point viewReflectedVector=operator-(dirVector,scalarProduct(normalVector,2.0*dotProduct(normalVector,dirVector)));
        viewReflectedVector=normalize(viewReflectedVector);

        phong+=scalingFactor*pow(dotProduct(viewReflectedVector,pointToLightVector),intersectedObject.shininess);



        
    }

    for(int i=0;i<spotLights.size();i++){
        Point pointToLightVector=operator-(spotLights[i].position,intersectedObject.intersectionPoint);
        pointToLightVector=normalize(pointToLightVector);
        
        Point lightToPointVector=operator-(intersectedObject.intersectionPoint,spotLights[i].position);
        lightToPointVector=normalize(lightToPointVector);

        Point normalVector=intersectedObject.normalVector;
        normalVector=normalize(normalVector);

        Point lightReflectedVector=operator-(lightToPointVector,scalarProduct(normalVector,2.0*dotProduct(normalVector,lightToPointVector)));

        Point sourceDirection=spotLights[i].direction;
        sourceDirection=normalize(sourceDirection);
        double dottemp = dotProduct(sourceDirection,lightToPointVector);
        double angle=acos(dottemp);
        angle=angle*180.0/M_PI;
        angle = 180-angle;
        if(angle>spotLights[i].cutoffAngle){
            
            continue;
        }

        double t1=(spotLights[i].position.x-intersectedObject.intersectionPoint.x)/pointToLightVector.x;

        IntersectedObject intersectedObject3=doesIlluminate(spotLights[i].position,lightToPointVector);

        double t2=intersectedObject3.t;

        if(((t2-t1)*(t2-t1))>0.0000001){
            continue;
        }

        double scalingFactor = exp(-t1*t1*spotLights[i].falloff);
        double NdotL=dotProduct(normalVector,pointToLightVector);
        NdotL < 0? NdotL=0:NdotL=NdotL;
        lambert+=scalingFactor*NdotL;
        
        Point viewReflectedVector=operator-(dirVector,scalarProduct(normalVector,2.0*dotProduct(normalVector,dirVector)));
        
        viewReflectedVector=normalize(viewReflectedVector);

        double LdotV=dotProduct(viewReflectedVector,pointToLightVector);
        LdotV < 0? LdotV=0:LdotV=LdotV;
        phong+=scalingFactor*pow(LdotV,intersectedObject.shininess);



        
    }

    


    Color ac=scalarProductColor(intersectedObject.color,intersectedObject.ka);
    Color dc=scalarProductColor(intersectedObject.color,intersectedObject.kd*lambert);
    Color sc=scalarProductColor(intersectedObject.color,intersectedObject.ks*phong);

    intersectedObject.color=operator+(ac,operator+(dc,sc));

    if(levelOfRecursion>0){
        Point viewReflectedVector=operator-(dirVector,scalarProduct(intersectedObject.normalVector,2.0*dotProduct(intersectedObject.normalVector,dirVector)));
        viewReflectedVector=normalize(viewReflectedVector);
        Point newIntersectionPoint=operator+(intersectedObject.intersectionPoint,scalarProduct(intersectedObject.normalVector,0.0001));
        IntersectedObject intersectedObject4=rayTrace(newIntersectionPoint,viewReflectedVector,levelOfRecursion-1);
        intersectedObject.color=operator+(intersectedObject.color,scalarProductColor(intersectedObject4.color,intersectedObject.kr));
    }





    return intersectedObject;
}

void capture(){
    windowHeight=(nearPlane*tan(fovY/2.0 *( M_PI/180.0)))*2.0;
    windowWidth=windowHeight*aspectRatio;
   
    double pixelWidth=windowWidth/numberOfPixels;
    double pixelHeight=windowHeight/numberOfPixels;
    double halfPixelWidth=pixelWidth/2.0;
    double halfPixelHeight=pixelHeight/2.0;
    double topLeftX=pos.x+l.x*nearPlane-r.x*windowWidth/2.0+u.x*windowHeight/2.0;
    double topLeftY=pos.y+l.y*nearPlane-r.y*windowWidth/2.0+u.y*windowHeight/2.0;
    double topLeftZ=pos.z+l.z*nearPlane-r.z*windowWidth/2.0+u.z*windowHeight/2.0;
    topLeftX=topLeftX+r.x*halfPixelWidth-u.x*halfPixelHeight;
    topLeftY=topLeftY+r.y*halfPixelWidth-u.y*halfPixelHeight;
    topLeftZ=topLeftZ+r.z*halfPixelWidth-u.z*halfPixelHeight;
    Point topLeftPoint(topLeftX,topLeftY,topLeftZ);
    int nearest;
    double t,tMin;
    
    bitmap_image image(numberOfPixels,numberOfPixels);
    for(int i=0;i<numberOfPixels;i++){
        for(int j=0;j<numberOfPixels;j++){

            



            double x=topLeftPoint.x+r.x*pixelWidth*j-u.x*pixelHeight*i;
            double y=topLeftPoint.y+r.y*pixelWidth*j-u.y*pixelHeight*i;
            double z=topLeftPoint.z+r.z*pixelWidth*j-u.z*pixelHeight*i;
            Point pixelPoint(x,y,z);

            double dirX=pixelPoint.x-pos.x;
            double dirY=pixelPoint.y-pos.y;
            double dirZ=pixelPoint.z-pos.z;
            Point dirVector(dirX,dirY,dirZ);
            dirVector=normalize(dirVector);

            
            IntersectedObject closest = rayTrace(pixelPoint, dirVector, levelOfRecursion);
            Color color = closest.color;
            
            image.set_pixel(j,i,(unsigned char)(color.r * 255.0),
                   (unsigned char)(color.g * 255.0),
                   (unsigned char)(color.b * 255.0));
            // pointBuffer[i][j][0]=color.r;
            // pointBuffer[i][j][1]=color.g;
            // pointBuffer[i][j][2]=color.b;


            
        }
        // printf("here 5\n");


    }
    image.save_image("out.bmp");

    cout <<"image captured" << endl;
            
   

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

        case '0':
            capture();
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
    int speed = 10;
	switch(key){
		case GLUT_KEY_UP:		//down arrow key
			pos.x+=l.x *speed;
			pos.y+=l.y *speed;
			pos.z+=l.z *speed;
			break;
		case GLUT_KEY_DOWN:		// up arrow key
			pos.x-=l.x*speed;
			pos.y-=l.y *speed;
			pos.z-=l.z*speed;
			break;

		case GLUT_KEY_RIGHT:
			pos.x+=r.x*speed;
			pos.y+=r.y*speed;
			pos.z+=r.z*speed;
			break;
		case GLUT_KEY_LEFT :
			pos.x-=r.x*speed;
			pos.y-=r.y*speed;
			pos.z-=r.z*speed;
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
    // gluPerspective(fovY, aspectRatio, nearPlane, farPlane);
    
    //glTranslatef(0.0f, 0.0f, -5.0f); // Move the octahedron back along the z-axis

    // control viewing (or camera)
    
    gluLookAt(pos.x,pos.y,pos.z,
              pos.x+l.x,pos.y+l.y,pos.z+l.z,
              u.x,u.y,u.z);

    // drawFirstCylinder();
    drawXYPlane(checkerboardWidth);
    drawAxis();
    glRotatef(rotateAngle, 0,1,0);
    for (Sphere& sphere : spheres) {
        sphere.drawSphere();
    }
    for(Cube& cube: cubes){
        cube.drawCube();
    }
    
    for(Pyramid& pyramid: pyramids){
        pyramid.drawPyramid();
    }

    for(pointLight& light: pointLights){
        light.drawPointLight();
    }

    for(spotLight& light: spotLights){
        

        light.drawSpotLight();
    }
    // drawAllCylinder();
    
    // glRotatef(-45,1,0,0);
// drawOctahedron();
//  drawAllSphere();
 
 

   


    glFlush();

}


/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
  
    // Set the viewport to cover the new window
     glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
  
    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(fovY, aspectRatio, nearPlane, farPlane);
}

int main(int argc, char** argv)
{
    pos.x=100;pos.y=0;pos.z=50;
    

    l.x=-1;l.y=0;l.z=0;
    
    u.x=0;u.y=0;u.z=1;
   
    r=crossProduct(l,u);
    


    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);
    glutInitWindowSize(600, 600);
    glutCreateWindow("Ray Tracing");
    glEnable(GL_DEPTH_TEST);
    readFromFile("description.txt");
    // buildUnitPositiveX(5);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);

    glutMainLoop();

    return 0;
}
