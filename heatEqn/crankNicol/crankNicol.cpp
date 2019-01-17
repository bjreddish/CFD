#include <iostream>
#include <fstream>
#include <math.h>

#define pi 3.14159265
// Declare variables
int numPoints=7;          //number of nodes
float timeEnd=0.75;        // end time
float timeCurr=0.0;        // current time - initialized at zero
float deltax;              // discritized length
float deltat=0.05;        // time step
float L=1;                 // length of our 1-D object 
float alpha=0.1;           // diffusivity constant


int main() 
{
    //***************************************//
    //******************SETUP****************//
    //***************************************//

	// Setting up output file
	std::ofstream outputFile;
	outputFile.open("output.txt");

	//Set up arrays
    //Our problem has numPoints number of nodes but out matrix will be
    // numPoints-1 x numPoints-1 matrix becuase that is how many unknows we have
	deltax = L/(numPoints-1);
    float A[numPoints-3];            //A = alpha*deltat/(2*deltax^2)   
    float B[numPoints-2];            //B = 1+ alpha*deltat/(2*deltax^2)
    float C[numPoints-3];             //C = A = alpha*deltat/(2*deltax^2)       
    float D[numPoints-2];                           // b values for out matrix of Ax=b
    float F[numPoints-2];                           // x vlaues for our matrix of Ax=b in this case we want T
    float C_star[numPoints-2];                      // new center diagonal 
    float D_star[numPoints-2];                      // new right diagonal

    // Init A
    for (int i=0;i<numPoints-3;i++)
    {
        A[i] = alpha*deltat/(2*deltax*deltax);
    }

    // Init B
    for (int i=0;i<numPoints-2;i++)
    {
        B[i] = 1+ alpha*deltat/(2*deltax*deltax);
    }

    // Init C
    for (int i=0;i<numPoints-3;i++)
    {
        C[i] = alpha*deltat/(2*deltax*deltax);
    }


    // Populate the initial conditions as a sin wave
    for (int i = 0; i <= numPoints-3; i++) 
    {
    	F[i] = sin(pi*(i+1)*deltax/L);
        //F = [node_1 ... node_n-1] 
        //node_0 node_n are bc's so they are not i the matrix
    }
    // Print out initial conditions for the sim t=0
    outputFile << "0";
    std::cout << "0";
    for (int i = 0; i <= numPoints-3; i++) 
    {
    	outputFile << "," << F[i];
        std::cout <<"," << F[i];
    }
    outputFile << ",0\n";
    std::cout <<  ",0\n";

    //***************************************//
    //******************SOLVE****************//
    //***************************************//

    // Itterate through the time steps
/*    while (timeCurr<timeEnd)
    {
    	outputFile << "0";// far left node
        // first 



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
    }*/

    outputFile.close();
    std::cout << "Done.";
    return 0;
}
