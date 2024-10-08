#include<iostream>
#include<GL/glut.h>
#include <cmath>
#include <ctime>
//glut windowr sathe related
//gl diye shuru func gula opengl er core function
//glu - gl diye shuru kichu utility function

double minAngle=0;
double secAngle=0;
double hourAngle=0;
double pendulumAngle=0;
int hours=0,minutes=0,seconds=0;
double theta=0.0;
double timeTracker=0.0; 
double omega=2*M_PI/2;
//int flag=0;
     

void angleInit(double sec, double min, double hour){
    secAngle=(360-sec*6+90)*M_PI/180;
    minAngle=(360-(min+sec/60)*6+90)*M_PI/180;
    hourAngle=(360-(hour+sec/3600+min/60)*30+90)*M_PI/180;
}

void middlePoint(){
    glPointSize(5);//point 10*10 pixels hobe, eta glBegin er aage deya lagbe, mustt
    //ekhan theke context shuru,ekhetre point draw er context, context bole kivabe draw, glEnd diye shesh
    glBegin(GL_POINTS);//() er majhe mode dite hobe, ki akte chai arki
        glColor3d(1,1,1);
        glVertex2d(0,0);//0,0 te point dekha jabe
    glEnd();
}

void secHandler(){
    //glLineSize(5);
    glColor3f(0.3,0.5,0.4);
    //  glBegin(GL_LINES);
    //     glVertex2d(0,0);
    //     glVertex2d(0.3*cos(secAngle),0.3*sin(secAngle));
    // glEnd();
    glBegin(GL_QUADS);
        glVertex2d(0,0);
        
        glVertex2d(0.25*cos(secAngle-0.03),0.25*sin(secAngle-0.03));
        glVertex2d(0.3*cos(secAngle),0.3*sin(secAngle));
        glVertex2d(0.25*cos(secAngle+0.03),0.25*sin(secAngle+0.03));
    

     glEnd();

}

void mintHandler(){
    //glLineSize(5);
    glColor3f(0.2,0.4,0.4);
    //  glBegin(GL_LINES);
    //     glVertex2d(0,0);
    //     glVertex2d(0.25*cos(minAngle),0.25*sin(minAngle));
    // glEnd();

        glBegin(GL_QUADS);
        glVertex2d(0,0);
        
        glVertex2d(0.20*cos(minAngle-0.04),0.20*sin(minAngle-0.04));
        glVertex2d(0.25*cos(minAngle),0.25*sin(minAngle));
        glVertex2d(0.20*cos(minAngle+0.04),0.20*sin(minAngle+0.04));
    

     glEnd();
}

void hourHandler(){
    //glLineSize(5);
     glColor3f(0.1,0.5,0.5);
    //  glBegin(GL_LINES);
    //     glVertex2d(0,0);
    //     glVertex2d(0.2*cos(hourAngle),0.2*sin(hourAngle));
    // glEnd();
      glBegin(GL_QUADS);
        glVertex2d(0,0);
        
        glVertex2d(0.12*cos(hourAngle-0.06),0.12*sin(hourAngle-0.06));
        glVertex2d(0.2*cos(hourAngle),0.2*sin(hourAngle));
        glVertex2d(0.12*cos(hourAngle+0.06),0.12*sin(hourAngle+0.06));
    

     glEnd();
}

void clockCircle(){
    float cx, cy, r,r1,r2;
     glLineWidth(1.5);
     //bairer circle 1
    glBegin(GL_LINE_LOOP);  // All vertices form a single loop of single pixel width
       
        glColor3f(0.5f,0.5f,0.8f);  // Light-blue
        cx = 0;
        cy = 0;
        r = 0.35;
        for (float theta = 0; theta < 360; theta += 10) {
            float x = cx + r * cos(theta/180*M_PI);
            float y = cy + r * sin(theta/180*M_PI);
            glVertex2f(x,y);
        }
    glEnd();
    //bairer circle 2
     glBegin(GL_LINE_LOOP);  // All vertices form a single loop of single pixel width
       
        glColor3f(0.5f,0.3f,0.4f);  // Light-blue
        cx = 0;
        cy = 0;
        r = 0.36;
        for (float theta = 0; theta < 360; theta += 10) {
            float x = cx + r * cos(theta/180*M_PI);
            float y = cy + r * sin(theta/180*M_PI);
            glVertex2f(x,y);
        }
    glEnd();
     //glColor3f(0.5f,0.5f,1.0f);
  
    //hour er vaag
    cx = 0;
    cy = 0;
    r1 = 0.35;
    r2 = 0.33;
    for (float theta = 0; theta < 360; theta += 6) {
            float x1 = cx + r1 * cos(theta/180*M_PI);
            float y1 = cy + r1 * sin(theta/180*M_PI);
            float x2 = cx + r2 * cos(theta/180*M_PI);
            float y2 = cy + r2 * sin(theta/180*M_PI);
            glBegin(GL_LINES);
            glColor3f(0.2f,0.5f,1.0f);
                glVertex2f(x1,y1);
                glVertex2f(x2,y2);
            glEnd();
        }
    //hour er vaag(3,6,9,12)
      cx = 0;
    cy = 0;
    r1 = 0.35;
    r2 = 0.31;
    for (float theta = 0; theta < 360; theta += 30) {
            float x1 = cx + r1 * cos(theta/180*M_PI);
            float y1 = cy + r1 * sin(theta/180*M_PI);
            float x2 = cx + r2 * cos(theta/180*M_PI);
            float y2 = cy + r2 * sin(theta/180*M_PI);
            glBegin(GL_LINES);
            glColor3f(0.5f,0.8f,1.0f);
                glVertex2f(x1,y1);
                glVertex2f(x2,y2);
            glEnd();
        }
}

void initPendulum(){
     theta=30.0*cos(omega*timeTracker);

    
     glColor3f(0.3,0.5,0.4);
    
    
    //  glBegin(GL_LINES);
    //  //glLineWidth(5);
    //     glVertex2d(0,-0.35);
    //     glVertex2d(0.75*cos((theta+270)*M_PI/180),0.75*sin((theta+270)*M_PI/180));
    // glEnd();

     glBegin(GL_QUADS);
        glVertex2d(0,-0.35);
        
        glVertex2d(0.5*cos((theta+270)*M_PI/180-0.05),0.5*sin((theta+270)*M_PI/180-0.05)-0.35);
        glVertex2d(0.52*cos((theta+270)*M_PI/180),0.52*sin((theta+270)*M_PI/180)-0.35);
        glVertex2d(0.50*cos((theta+270)*M_PI/180+0.05),0.50*sin((theta+270)*M_PI/180+0.05)-0.35);
    

     glEnd();
    float cx, cy, r;
     //glLineWidth(2);
    glBegin(GL_POLYGON);  // All vertices form a single loop of single pixel width
       
        glColor3f(0.4f,1.0f,1.0f);  // Light-blue
        cx = 0.55*cos((theta+270)*M_PI/180);
        cy = 0.55*sin((theta+270)*M_PI/180)-0.35;
        r = 0.03;
        for (float theta = 0; theta < 360; theta += 10) {
            float x = cx + r * cos(theta/180*M_PI);
            float y = cy + r * sin(theta/180*M_PI);
            glVertex2f(x,y);
        }
    glEnd();

}

void timer(int value) {
    glutPostRedisplay();    // Post re-paint request to activate display()
    glutTimerFunc(10, timer, 0); // Call next 'timer' milliseconds later
}

void update(int value) {
    

    ///write angle change code
    secAngle=secAngle-6*M_PI/180;
    minAngle=minAngle-M_PI*0.1/180;
    hourAngle=hourAngle-M_PI/(120*180);
    //printf("%lf\n",secAngle);
    glutPostRedisplay(); 
        glutTimerFunc(1000, update, 0); // Call next 'timer' milliseconds later
}



void pendulumUpdate(int value){

    
     timeTracker+=0.02;
    glutPostRedisplay(); 
    glutTimerFunc(20, pendulumUpdate, 0);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);//color clear kore
   
    secHandler();
    mintHandler();
    hourHandler();
    clockCircle();
    middlePoint();
    initPendulum();
    //pointer vertex
   
    glFlush();//render now, etokkhon ja likhsi ta buffer theke ene dekhabe
}//display te ki dekhte chai ta thakbe

int main(int argc, char** argv){
    glutInit(&argc,argv);
    glutInitWindowSize(640,640);//window resolution arki
    glutInitWindowPosition(100,50);//initially window koi thakbe
    glutCreateWindow("Clock");
    // Get the current time
    std::time_t currentTime = std::time(nullptr);

    // Convert the current time to the local time
    std::tm* localTime = std::localtime(&currentTime);

    // Extract hours, minutes, and seconds
    hours = localTime->tm_hour;
    minutes = localTime->tm_min;
    seconds = localTime->tm_sec;
    angleInit(seconds,minutes,hours);
    glutDisplayFunc(display);//display render korte call, like display window minimize kore then call korlam, tokhon call hobe eta
    // glutMainLoop();//infifnite loop, event queue theke next event execute korbe, eta ekebare last e add
    
    glutTimerFunc(0, timer, 0);             // First timer call immediately
    glutTimerFunc(0, update, 0);
    glutTimerFunc(0,pendulumUpdate,0);
    glutMainLoop(); 
    return 0;

}