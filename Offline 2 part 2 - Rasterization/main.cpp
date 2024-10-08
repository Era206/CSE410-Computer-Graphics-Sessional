#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stack>
#include <iomanip> 
#include <cmath>
#include "bitmap_image.hpp"
using namespace std;


double eyeX, eyeY, eyeZ;
double lookX, lookY, lookZ;
double upX, upY, upZ;
double fovY, aspectRatio, nearPlane, farPlane;

 vector< vector<double>> M(4, vector<double>(4, 0.0));
 stack< vector< vector<double>>> matrixStack;
 stack< vector< vector<double>>> stage1MatrixStack;
 stack< vector< vector<double>>> stage2MatrixStack;
 stack< vector<double>> stage1PointStack;



struct Point {
    double x;
    double y;
    double z;

    Point(double xVal, double yVal, double zVal) : x(xVal), y(yVal), z(zVal) {}
};

struct Color{
    double r;
    double g;
    double b;

    Color(double rVal, double gVal, double bVal) : r(rVal), g(gVal), b(bVal) {}
};

vector<Color> stage3ColorVector;

vector< vector<Point>> stage3TriangleVector;

static unsigned long int g_seed = 1; 
inline int randomValueAssign() 
{ 
    g_seed = (214013 * g_seed + 2531011); 
    return (g_seed >> 16) & 0x7FFF; 
} 

// Function to calculate the cross product of two 3D position vectors
Point crossProduct(const Point& v1, const Point& v2) {
    Point result(0.0,0.0,0.0);
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}

// Function to calculate the dot product of two 3D position vectors
double dotProduct(const Point& v1, const Point& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// Function to calculate the vector R(x, a, angle)
Point R(const Point& x, const Point& a, double angle) {
    angle = angle * M_PI / 180.0;

    double cosAngle =  cos(angle);
    double sinAngle =  sin(angle);

    double dotProductResult = dotProduct(a, x);

    Point crossProductResult = crossProduct(a, x);

    Point result(0.0,0.0,0.0);
    result.x = cosAngle * x.x + (1 - cosAngle) * dotProductResult * a.x + sinAngle * crossProductResult.x;
    result.y = cosAngle * x.y + (1 - cosAngle) * dotProductResult * a.y + sinAngle * crossProductResult.y;
    result.z = cosAngle * x.z + (1 - cosAngle) * dotProductResult * a.z + sinAngle * crossProductResult.z;

    return result;
}


// Function to perform matrix multiplication of two 2D vector matrices
 vector< vector<double>> matrixMultiplication(const  vector< vector<double>>& matrix1,const  vector< vector<double>>& matrix2) {
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int rows2 = matrix2.size();
    int cols2 = matrix2[0].size();

    if (cols1 != rows2) {
         cout << "Matrix dimensions are not compatible for multiplication." <<  endl;
        return  vector< vector<double>>(); 
    }

     vector< vector<double>> result(rows1,  vector<double>(cols2, 0.0));

    // Perform matrix multiplication
    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols2; ++j) {
            for (int k = 0; k < cols1; ++k) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return result;
}

// Function to write a specified number of rows of the matrix to a text file
void writeMatrixToFile(const  vector< vector<double>>& matrix, int numRows, const  string& filename) {
      ofstream outputFile(filename,  ios::app); // Open the file in append mode
    
    if (!outputFile.is_open()) {
         cerr << "Error opening file: " << filename <<  endl;
        return;
    }

    
    outputFile <<  fixed <<  setprecision(7);

    
    outputFile << matrix[0][0] << " " << matrix[1][0] << " " << matrix[2][0];
        outputFile <<  endl;

    //outputFile.close();
}





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

void stage1Func(const  string& filename) {
     ifstream inputFile(filename);
    if (!inputFile.is_open()) {
         cout << "Error opening file: " << filename <<  endl;
        return;
    }

    // Read fixed parameters
    for (int i = 0; i < 4; ++i) {
         string line;
         getline(inputFile, line);
         vector<double> values = parseLine(line);
        if(i==0){
            eyeX=values[0];
            eyeY=values[1];
            eyeZ=values[2];
        }
        if(i==1){
            lookX=values[0];
            lookY=values[1];
            lookZ=values[2];
        }
        if(i==2){
            upX=values[0];
            upY=values[1];
            upZ=values[2];
        }
        if(i==3){
            fovY=values[0], aspectRatio=values[1], nearPlane=values[2], farPlane=values[3];
        }
    }

     string line;
    while ( getline(inputFile, line)) {
        if (line == "end") {
            break; 
        }

        if (line == "triangle") {
            
          
            for (int i = 0; i < 3; ++i) {
                //printf("i: %d\n",i);
                 getline(inputFile, line);
                 vector<double> values = parseLine(line);
                 vector< vector<double>> triangleMatrix(4,  vector<double>(1, 0.0));
                 vector< vector<double>> resultMatrix(4,  vector<double>(1, 0.0));

                for (int j = 0; j < 3; ++j) {
                    triangleMatrix[j][0] = values[j];
                }
                triangleMatrix[3][0] = 1.0;

                resultMatrix=matrixMultiplication(M,triangleMatrix);
                int numRowsToWrite = resultMatrix.size() - 1;

                
                writeMatrixToFile(resultMatrix, numRowsToWrite, "stage1.txt");
                
                stage1MatrixStack.push(resultMatrix);
              

            }

             ofstream outputFile("stage1.txt",  ios::app);
            outputFile <<  endl;
           
            
        } else if (line == "translate") {
            // Read tx, ty, tz parameters (1 line)
             getline(inputFile, line);
             vector<double> values = parseLine(line);
             vector< vector<double>> translateMatrix(4,  vector<double>(4, 0.0));
            
            for (int j = 0; j < 3; ++j) {
                    translateMatrix[j][3] = values[j];
                    translateMatrix[j][j] = 1.0;
                }
                translateMatrix[3][3] = 1.0;

                M=matrixMultiplication(M,translateMatrix);

                // // Accessing and printing the matrix
                // for (int i = 0; i < 4; ++i) {
                //     for (int j = 0; j < 4; ++j) {
                //          cout << M[i][j] << " ";
                //     }
                //      cout <<  endl;
                // }
            
        } else if (line == "scale") {
            // Read sx, sy, sz parameters (1 line)
             getline(inputFile, line);
             vector<double> values = parseLine(line);
            // Assign values to the matrix (using scaling matrix)
             vector< vector<double>> scaleMatrix(4,  vector<double>(4, 0.0));
            
            for (int j = 0; j < 3; ++j) {
                    scaleMatrix[j][j] = values[j];
                }
                scaleMatrix[3][3] = 1.0;

                M=matrixMultiplication(M,scaleMatrix);


                // // Accessing and printing the matrix
                // for (int i = 0; i < 4; ++i) {
                //     for (int j = 0; j < 4; ++j) {
                //          cout << M[i][j] << " ";
                //     }
                //      cout <<  endl;
                // }
            // ...
        } else if (line == "rotate") {
            // Read angle, ax, ay, az parameters (1 line)
             getline(inputFile, line);
             vector<double> values = parseLine(line);
            double ax, ay, az;
            ax=values[1];
            ay=values[2];
            az=values[3];

            Point a(ax, ay, az);
            double magA =  sqrt(ax * ax + ay * ay + az * az);
            a.x /= magA;
            a.y /= magA;
            a.z /= magA;


            double angle=values[0];

            // Calculate c1, c2, c3
            Point c1 = R(Point(1.0, 0.0, 0.0), a, angle);
            Point c2 = R(Point(0.0, 1.0, 0.0), a, angle);
            Point c3 = R(Point(0.0, 0.0, 1.0), a, angle);
            // Assign values to the matrix (using rotation matrix)
             vector< vector<double>> rotateMatrix(4,  vector<double>(4, 0.0));
            
            rotateMatrix[0][0] = c1.x;
            rotateMatrix[0][1] = c2.x;
            rotateMatrix[0][2] = c3.x;
            rotateMatrix[1][0] = c1.y;
            rotateMatrix[1][1] = c2.y;
            rotateMatrix[1][2] = c3.y;
            rotateMatrix[2][0] = c1.z;
            rotateMatrix[2][1] = c2.z;
            rotateMatrix[2][2] = c3.z;
            

            rotateMatrix[3][3] = 1.0;

            M=matrixMultiplication(M,rotateMatrix);
            
        } else if (line == "push") {
            

            matrixStack.push(M);
        } else if (line == "pop") {
            

            M=matrixStack.top();
            matrixStack.pop();
        } else {
            
             cout << "Unknown command: " << line <<  endl;
        }
    }

    inputFile.close();
}


void stage2Func(){
    Point look(lookX, lookY, lookZ);
    Point eye(eyeX, eyeY, eyeZ);
    Point l(0.0,0.0,0.0);
    l.x=look.x-eye.x;
    l.y=look.y-eye.y;
    l.z=look.z-eye.z;

    double magL =  sqrt(l.x * l.x + l.y * l.y + l.z * l.z);
    l.x /= magL;   
    l.y /= magL;
    l.z /= magL;

    Point up(upX, upY, upZ);
    Point r=crossProduct(l,up);
    double magR =  sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
    r.x /= magR;
    r.y /= magR;
    r.z /= magR;

    Point u=crossProduct(r,l);
     vector< vector<double>> T(4,  vector<double>(4, 0.0));
    for(int i=0;i<4;i++){
        T[i][i]=1.0;
    }
    T[0][3]=-eye.x;
    T[1][3]=-eye.y;
    T[2][3]=-eye.z;

     vector< vector<double>> R(4,  vector<double>(4, 0.0));
    R[0][0]=r.x;
    R[0][1]=r.y;
    R[0][2]=r.z;
    R[1][0]=u.x;
    R[1][1]=u.y;
    R[1][2]=u.z;
    R[2][0]=-l.x;
    R[2][1]=-l.y;
    R[2][2]=-l.z;
    R[3][3]=1.0;

     vector< vector<double>> V(4,  vector<double>(4, 0.0));
    V=matrixMultiplication(R,T);
    int size=stage1MatrixStack.size();
     stack< vector< vector<double>>> tempVector;
   // Reversing the stack of vectors
    while (!stage1MatrixStack.empty()) {
        tempVector.push(stage1MatrixStack.top());
        stage1MatrixStack.pop();
    }

    stage1MatrixStack =  move(tempVector);

    for(int i=0;i<size;i++){
         vector< vector<double>> triangleMatrix(4,  vector<double>(1, 0.0));
         vector< vector<double>> resultMatrix(4,  vector<double>(1, 0.0));

        triangleMatrix=stage1MatrixStack.top();
        stage1MatrixStack.pop();
        triangleMatrix[3][0] = 1.0;
        resultMatrix=matrixMultiplication(V,triangleMatrix);
        int numRowsToWrite = resultMatrix.size() - 1;

    
        writeMatrixToFile(resultMatrix, numRowsToWrite, "stage2.txt");
        stage2MatrixStack.push(resultMatrix);

        if((i+1)%3==0){
                ofstream outputFile("stage2.txt",  ios::app);
            outputFile <<  endl;
        }
    }

}


void stage3Func(){
    double fovX=fovY*aspectRatio;
    double t=nearPlane* tan(fovY/2.0*M_PI/180.0);
    double r=nearPlane* tan(fovX/2.0*M_PI/180.0);

     vector< vector<double>> P(4,  vector<double>(4, 0.0));
    P[0][0]=nearPlane/r;
    P[1][1]=nearPlane/t;
    P[2][2]=-(farPlane+nearPlane)/(farPlane-nearPlane);
    P[2][3]=-(2*farPlane*nearPlane)/(farPlane-nearPlane);
    P[3][2]=-1.0;

    

    int size=stage2MatrixStack.size();
     stack< vector< vector<double>>> tempVector;
  
    while (!stage2MatrixStack.empty()) {
        tempVector.push(stage2MatrixStack.top());
        
        stage2MatrixStack.pop();
    }

    Point p(0.0,0.0,0.0);
    vector<Point> tempTriangleVector;

    for(int i=0;i<size;i++){
         vector< vector<double>> triangleMatrix(4,  vector<double>(1, 0.0));
                 vector< vector<double>> resultMatrix(4,  vector<double>(1, 0.0));


                triangleMatrix=tempVector.top();
                tempVector.pop();
               

        resultMatrix=matrixMultiplication(P,triangleMatrix);
                int numRowsToWrite = resultMatrix.size() - 1;
                


                // if(resultMatrix[3][0]!=1.0){
                    resultMatrix[0][0]=resultMatrix[0][0]/resultMatrix[3][0];
                    resultMatrix[1][0]=resultMatrix[1][0]/resultMatrix[3][0];
                    resultMatrix[2][0]=resultMatrix[2][0]/resultMatrix[3][0];
                    resultMatrix[3][0]=1.0;

                // }
                
                p.x=resultMatrix[0][0];
                p.y=resultMatrix[1][0];
                p.z=resultMatrix[2][0];
                tempTriangleVector.push_back(p);
                
                writeMatrixToFile(resultMatrix, numRowsToWrite, "stage3.txt");
                

                if((i+1)%3==0){
                    int r = randomValueAssign();
                    int g = randomValueAssign();
                    int b = randomValueAssign();
                    Color color(r, g, b);
                    stage3ColorVector.push_back(color);
                    vector<Point> tempVector2;
                    tempVector2=move(tempTriangleVector);
                    stage3TriangleVector.push_back(tempVector2);
                    tempTriangleVector.clear();
                     ofstream outputFile("stage3.txt",  ios::app);
                    outputFile <<  endl;
                }
                
    }

    // int i=0;
    // for (const auto& vector : stage3TriangleVector) {
    //     cout << "Vector of PosVector:" << endl;
    //     for (const Point& pos : vector) {
    //         cout << "  PosVector: (" << pos.x << ", " << pos.y << ", " << pos.z << ")" << endl;
    //     }
    //     cout <<"color of posvector:"<<stage3ColorVector[i].r<<" "<<stage3ColorVector[i].g<<" "<<stage3ColorVector[i].b<<endl;
    //     i++;
    // }



}




struct Plane{
    double a;
    double b;
    double c;
    double d;

    Plane(double aVal, double bVal, double cVal, double dVal) : a(aVal), b(bVal), c(cVal), d(dVal) {}
};


Plane planeEquation(const Point& p1, const Point& p2, const Point& p3) {
    Point v1(p2.x-p1.x, p2.y-p1.y, p2.z-p1.z);
    Point v2(p3.x-p1.x, p3.y-p1.y, p3.z-p1.z);
    Point normalVector=crossProduct(v1,v2);
    double d=-(normalVector.x*p1.x+normalVector.y*p1.y+normalVector.z*p1.z);
    Plane result(normalVector.x, normalVector.y, normalVector.z, d);
    return result;
}


double zValue(const Plane& plane, const Point& p){
    double z=(-plane.a*p.x-plane.b*p.y-plane.d)/plane.c;
    return z;
}


double*** zBufferWithColor;


void initializeZBufferWithColor(int dim1, int dim2) {
    zBufferWithColor = new double**[dim1];
    for (int i = 0; i < dim1; ++i) {
        zBufferWithColor[i] = new double*[dim2];
        for (int j = 0; j < dim2; ++j) {
            zBufferWithColor[i][j] = new double[4];
            zBufferWithColor[i][j][0] = 1.0;  // zVal
            zBufferWithColor[i][j][1] = 0.0;  // r
            zBufferWithColor[i][j][2] = 0.0;  // g
            zBufferWithColor[i][j][3] = 0.0;  // b
        }
    }
}


void deleteZBufferWithColor(int dim1, int dim2) {
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            delete[] zBufferWithColor[i][j];
        }
        delete[] zBufferWithColor[i];
    }
    delete[] zBufferWithColor;
}


double min(double a, double b, double c){
    double min=a;
    if(b<min){
        min=b;
    }
    if(c<min){
        min=c;
    }
    return min;
}


double max(double a, double b, double c){
    double max=a;
    if(b>max){
        max=b;
    }
    if(c>max){
        max=c;
    }
    return max;
}

int screenWidth, screenHeight;
double rightLimitX=1.0, leftLimitX=-1.0, topLimitY=1.0, bottomLimitY=-1.0;
double posLimitZ=1.0, negLimitZ=-1.0;
double dx,dy;
double TopY, LeftX;

void stage4Func(const  string& filename){
     ifstream inputFile(filename);
    if (!inputFile.is_open()) {
         cout << "Error opening file: " << filename <<  endl;
        return;
    }

   
    string line;
    getline(inputFile, line);
    vector<double> values = parseLine(line);
    
    screenWidth=values[0];
    screenHeight=values[1];

    initializeZBufferWithColor(screenHeight, screenWidth);
    dx=(rightLimitX-leftLimitX)/(screenWidth*1.0);
    dy=(topLimitY-bottomLimitY)/(screenHeight*1.0);
    TopY=topLimitY-dy/2.0;
    LeftX=leftLimitX+dx/2.0;
    int length=stage3TriangleVector.size();
    for(int i=0;i<length;i++){
        vector<Point> tempVector;
        tempVector=stage3TriangleVector[i];
        Point p1=tempVector[0];
        Point p2=tempVector[1];
        Point p3=tempVector[2];
        
        Plane plane=planeEquation(p1,p2,p3);
        double minY=min(p1.y,p2.y,p3.y);
        double maxY=max(p1.y,p2.y,p3.y);
        int rowMin=ceil((TopY-maxY)/dy);
        int rowMax=floor((TopY-minY)/dy);
        if(rowMin<0){
            rowMin=0;
        }
        if(rowMax>=screenHeight){
            rowMax=screenHeight-1;
        }
        for(int j=rowMin;j<=rowMax;j++){
            double y=TopY-j*dy;
            Point tempP1(0.0,0.0,0.0);
            Point tempP2(0.0,0.0,0.0);
            Point tempP3(0.0,0.0,0.0);
            if(p1.y==p2.y||p2.y==p3.y||p1.y==p3.y){
                p1.y==p2.y?(tempP1=p3,tempP2=p1,tempP3=p2):(p2.y==p3.y?(tempP1=p1,tempP2=p2,tempP3=p3):(tempP1=p2,tempP2=p3,tempP3=p1));
            }
            else{
                double ty12=(y-p1.y)/(p2.y-p1.y);
                double ty23=(y-p2.y)/(p3.y-p2.y);
                double ty31=(y-p3.y)/(p1.y-p3.y);
                if(ty12>=0&&ty12<=1 && ty23>=0&&ty23<=1){
                    tempP1=p2;
                    tempP2=p3;
                    tempP3=p1;
                }
                else if(ty12>=0&&ty12<=1 && ty31>=0&&ty31<=1){
                    tempP1=p1;
                    tempP2=p2;
                    tempP3=p3;
                }
                else if(ty23>=0&&ty23<=1 && ty31>=0&&ty31<=1){
                    tempP1=p3;
                    tempP2=p1;
                    tempP3=p2;
                }
            }
            

            double tTempP12=(y-tempP1.y)/(tempP2.y-tempP1.y);
            double tTempP13=(y-tempP1.y)/(tempP3.y-tempP1.y);
            double x1=tempP1.x+tTempP12*(tempP2.x-tempP1.x);
            double x2=tempP1.x+tTempP13*(tempP3.x-tempP1.x);
            double colMin=round((x1-LeftX)/dx);
            double colMax=round((x2-LeftX)/dx);
            

            if(colMin<0){
                colMin=0;
            }
            if(colMax>=screenWidth){
                colMax=screenWidth-1;
            }
            for(int k=colMin;k<=colMax;k++){
                double x=LeftX+k*dx;
                Point p(x,y,0.0);
                double z=zValue(plane,p);
                
                if(z>=negLimitZ && z<=posLimitZ){
                    if(z<zBufferWithColor[j][k][0]){
                        zBufferWithColor[j][k][0]=z;
                        zBufferWithColor[j][k][1]=stage3ColorVector[i].r;
                        zBufferWithColor[j][k][2]=stage3ColorVector[i].g;
                        zBufferWithColor[j][k][3]=stage3ColorVector[i].b;
                    }
                }

            }
            
        }
    }

 
    ofstream outputFile("z_buffer.txt",  ios::trunc);
    if (!outputFile.is_open()) {
         cout << "Error opening file: " << "z_buffer.txt" <<  endl;
        return;
    }
    outputFile <<  fixed <<  setprecision(6);
    
    for(int i=0;i<screenHeight;i++){
        for(int j=0;j<screenWidth;j++){
            if(zBufferWithColor[i][j][0]<1.0){
                outputFile<<zBufferWithColor[i][j][0]<<"\t";
            }
            
        }
        outputFile<<endl;
    }

    
    outputFile.close();
    bitmap_image image(screenWidth, screenHeight);
    
    for(int i=0;i<screenHeight;i++){
        for(int j=0;j<screenWidth;j++){
            // if(zBufferWithColor[i][j][0]<1.0){
                image.set_pixel(j,i,zBufferWithColor[i][j][1],zBufferWithColor[i][j][2],zBufferWithColor[i][j][3]);
            // }
            
        }
    }
    
    image.save_image("out.bmp");
    image.clear();



    

    
}






int main(int argc, char** argv) {
     ofstream outputFile("stage1.txt",  ios::trunc);
    outputFile.close();
     ofstream outputFile2("stage2.txt",  ios::trunc);
    outputFile2.close();
     ofstream outputFile3("stage3.txt",  ios::trunc);
    outputFile3.close();



    for (int i = 0; i < 4; ++i) {
        M[i][i] = 1.0;
    }

    
     string filename = "scene.txt";
     string filename2 = "config.txt";

    // Call the readFromFile function and pass the filename as an argument
    stage1Func(filename);
    stage2Func();
    stage3Func();
    stage4Func(filename2);
    deleteZBufferWithColor(screenHeight, screenWidth);
    

    return 0;
}

