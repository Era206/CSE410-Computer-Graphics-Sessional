#include <GL/glut.h>

const int gridSize = 10;  // Number of grid lines in each direction

void drawGrid() {
    glColor3f(0.5f, 0.5f, 0.5f);  // Set color to gray

    // Draw horizontal lines
    for (int i = -gridSize; i <= gridSize; ++i) {
        glBegin(GL_LINES);
        glVertex2f(-gridSize, i);
        glVertex2f(gridSize, i);
        glEnd();
    }

    // Draw vertical lines
    for (int i = -gridSize; i <= gridSize; ++i) {
        glBegin(GL_LINES);
        glVertex2f(i, -gridSize);
        glVertex2f(i, gridSize);
        glEnd();
    }
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluOrtho2D(-gridSize, gridSize, -gridSize, gridSize);

    drawGrid();

    glFlush();
    glutSwapBuffers();
}

void reshape(int width, int height) {
    glViewport(0, 0, width, height);
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Grid Drawing");

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);  // Set clear color to white

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);

    glutMainLoop();
    return 0;
}
