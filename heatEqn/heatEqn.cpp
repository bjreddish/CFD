#include <iostream>
#include <fstream>
#include <math.h>

#define pi 3.14159265
// Declare variables
int numPoints=20;  //number of elements to calculate
float timeEnd=0.75;    // end time
float time=0;
float bcLeft=0; // temperature bc on the left
float bcRight=0; // temperature bc on the right
float deltax;
float deltat=0.005;
float L=1;
float alpha=0.1; // diffusivity constant
float r;
float r2;

int main() 
{
	std::cout << "Beginning heat calculations\n";
	// Setting up output file
	std::ofstream outputFile;
	outputFile.open("output.txt");
	//Set up arrays
	deltax = L/numPoints;
    float tempOld[numPoints];
    float tempNew[numPoints];

    tempOld[0] = bcLeft;
    tempOld[numPoints-1] = bcRight;
    // Populate the initial conditions
    for (int j = 1; j <= numPoints-2; j++) // itterating from index 1 to 98 (index 0 and 99 are const)
    {
    	tempOld[j] = sin(pi*j*deltax/L);
    }
    // Print out initial conditions for the sim t=0
    outputFile << "0";
    for (int i = 1; i < numPoints; i++) 
    {
    	outputFile << "," << tempOld[i];
    }
    outputFile << "\n";

    r = alpha*deltat/(deltax*deltax);
    r2 = 1 - 2*r;
    while (time<timeEnd)
    {
    	outputFile << "0";// far left node
    	for (int j = 1; j <= numPoints-2; j++) // itterating from index 1 to 98 (index 0 and 99 are const)
   		{
    		tempNew[j] = tempOld[j] + alpha*deltat/(deltax*deltax)*(tempOld[j+1]-2*tempOld[j]+tempOld[j-1]);
    		outputFile << "," << tempNew[j];
    	}
    	outputFile << "," <<"0\n";// far right node
    	
    	for (int k =1; k<= numPoints-2;k++)
    	{
    		tempOld[k]=tempNew[k];
    	}	
    	time = time+deltat;
    }

    outputFile.close();
    std::cout << "Done.";
    return 0;
}