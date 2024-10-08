#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stack>
#include <iomanip> 
#include <cmath>
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



struct PosVector {
    double x;
    double y;
    double z;

    PosVector(double xVal, double yVal, double zVal) : x(xVal), y(yVal), z(zVal) {}
};

// Function to calculate the cross product of two 3D position vectors
PosVector crossProduct(const PosVector& v1, const PosVector& v2) {
    PosVector result(0.0,0.0,0.0);
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}

// Function to calculate the dot product of two 3D position vectors
double dotProduct(const PosVector& v1, const PosVector& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// Function to calculate the vector R(x, a, angle)
PosVector R(const PosVector& x, const PosVector& a, double angle) {
    angle = angle * M_PI / 180.0;

    double cosAngle =  cos(angle);
    double sinAngle =  sin(angle);

    double dotProductResult = dotProduct(a, x);

    PosVector crossProductResult = crossProduct(a, x);

    PosVector result(0.0,0.0,0.0);
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
            // Assign values to the matrix (using translation matrix)
            // ...
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

            PosVector a(ax, ay, az);
            double magA =  sqrt(ax * ax + ay * ay + az * az);
            a.x /= magA;
            a.y /= magA;
            a.z /= magA;


            double angle=values[0];

            // Calculate c1, c2, c3
            PosVector c1 = R(PosVector(1.0, 0.0, 0.0), a, angle);
            PosVector c2 = R(PosVector(0.0, 1.0, 0.0), a, angle);
            PosVector c3 = R(PosVector(0.0, 0.0, 1.0), a, angle);
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
            // ...
        } else if (line == "push") {
            // Implement push matrix operation (if you need to handle a stack of matrices)

            matrixStack.push(M);
        } else if (line == "pop") {
            // Implement pop matrix operation (if you need to handle a stack of matrices)

            M=matrixStack.top();
            matrixStack.pop();
        } else {
            // Unknown command or invalid format, you can handle this case accordingly
             cout << "Unknown command: " << line <<  endl;
        }
    }

    inputFile.close();
}


void stage2Func(){
    PosVector look(lookX, lookY, lookZ);
    PosVector eye(eyeX, eyeY, eyeZ);
    PosVector l(0.0,0.0,0.0);
    l.x=look.x-eye.x;
    l.y=look.y-eye.y;
    l.z=look.z-eye.z;

    double magL =  sqrt(l.x * l.x + l.y * l.y + l.z * l.z);
    l.x /= magL;   
    l.y /= magL;
    l.z /= magL;

    PosVector up(upX, upY, upZ);
    PosVector r=crossProduct(l,up);
    double magR =  sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
    r.x /= magR;
    r.y /= magR;
    r.z /= magR;

    PosVector u=crossProduct(r,l);
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
                
                writeMatrixToFile(resultMatrix, numRowsToWrite, "stage3.txt");
                

                if((i+1)%3==0){
                     ofstream outputFile("stage3.txt",  ios::app);
                    outputFile <<  endl;
                }
    }



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

    // Call the readFromFile function and pass the filename as an argument
    stage1Func(filename);
    stage2Func();
    stage3Func();

    return 0;
}

