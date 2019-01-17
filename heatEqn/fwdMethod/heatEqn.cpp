#include <iostream>
#include <fstream>
#include <math.h>

#define pi 3.14159265
// Declare variables
int numPoints=20;   //number of elements to calculate
float timeEnd=0.75; // end time
float timeCurr=0.0; // current time - initialized at zero
float bcLeft=0;     // temperature bc on the left
float bcRight=0;    // temperature bc on the right
float deltax;       // discritized length
float deltat=0.005; // time step
float L=1;          // length of our 1-D object 
float alpha=0.1;    // diffusivity constant

int main() 
{
	// Setting up output file
	std::ofstream outputFile;
	outputFile.open("output.txt");

	//Set up arrays
	deltax = L/numPoints;
    float tempOld[numPoints];
    float tempNew[numPoints];
    // Setting BC's on far left and far right
    tempOld[0] = bcLeft;
    tempOld[numPoints-1] = bcRight;

    // Populate the initial conditions as a sin wave
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


    // Itterate through the time steps
    while (timeCurr<timeEnd)
    {
    	outputFile << "0";// far left node
    	for (int j = 1; j <= numPoints-2; j++) // itterating from index 1 to 98 (index 0 and 99 are const)
   		{
    		tempNew[j] = tempOld[j] + alpha*deltat/(deltax*deltax)*(tempOld[j+1]-2*tempOld[j]+tempOld[j-1]);
    		outputFile << "," << tempNew[j];
    	}
    	outputFile << "," <<"0\n";// far right node
    	// Replace old array with new array for next step of the itteration
    	for (int k =1; k<= numPoints-2;k++)
    	{
    		tempOld[k]=tempNew[k];
    	}	
    	timeCurr = timeCurr+deltat;
    }

    outputFile.close();
    std::cout << "Done.";
    return 0;
}
