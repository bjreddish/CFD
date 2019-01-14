#include <iostream>
#include <fstream>
#include <math.h>

#define pi 3.14159265
// Declare variables
int numPoints=100;  //number of elements to calculate
int timeEnd=10;    // end time
float bcLeft=0; // temperature bc on the left
float bcRight=0; // temperature bc on the right
float deltax;
float deltat=0.01;
float L=1;


int main() 
{
	std::cout << "Beginning heat calculations\n";
	// Setting up output file
	std::ofstream outputFile;
	outputFile.open("output.txt");
	outputFile << 1;
	for (int i = 2; i <= numPoints; i++) 
    {
    	outputFile << "," <<i;
    }
	//Set up arrays
	deltax = L/numPoints;
    float tempOld[numPoints];
    float tempNew[numPoints];
    float temp[numPoints];
    temp[0] = bcLeft;
    temp[numPoints-1] = bcRight;
    // Populate the initial conditions
    for (int j = 1; j <= numPoints-2; j++) // itterating from index 1 to 98 (index 0 and 99 are const)
    {
    	temp[j] = sin(pi*j*deltax/L);
    }
    
    for (int i = 0; i < numPoints; i++) 
    {
    	std::cout << i << " : "<<temp[i] << "\n";
    }

    outputFile.close();
    std::cout << "Done.";
    return 0;
}