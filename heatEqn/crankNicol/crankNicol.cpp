#include <iostream>
#include <fstream>
#include <math.h>

#define pi 3.14159265
// Declare variables
int numPoints=7;          //number of nodes
float timeEnd=0.75;        // end time
float timeCurr=0.0;        // current time - initialized at zero
float deltax;              // discritized length
float deltat=0.001;        // time step
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

    std::cout << "Number of nodes: " << numPoints << "\n";
    std::cout << "Simulation time: " << timeEnd << " seconds\n";
    std::cout << "Time step: " << deltat << " seconds\n";

	//Set up arrays
    //Our problem has numPoints number of nodes but out matrix will be
    // numPoints-1 x numPoints-1 matrix becuase that is how many unknows we have
	deltax = L/(numPoints-1);
    std::cout << "Delta x: " << deltax << "\n";
    float A[numPoints-3];            //A = alpha*deltat/(2*deltax^2)   
    float B[numPoints-2];            //B = 1+ alpha*deltat/(2*deltax^2)
    float C[numPoints-3];             //C = A = alpha*deltat/(2*deltax^2)       
    float D[numPoints-2];                           // b values for out matrix of Ax=b
    float F[numPoints];                           // x vlaues for our matrix of Ax=b in this case we want T
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
        B[i] = -(1+ alpha*deltat/(2*deltax*deltax));
    }

    // Init C
    for (int i=0;i<numPoints-3;i++)
    {
        C[i] = alpha*deltat/(2*deltax*deltax);
    }


    // Populate the initial conditions as a sin wave
    for (int i = 0; i < numPoints; i++) 
    {
    	F[i] = sin(pi*(i)*deltax/L);
    }
    F[numPoints-1]=0;

    std::cout << "Size of A coefficient array: " << (sizeof(A)/sizeof(*A)) << "\n";
    std::cout << "Size of B coefficient array: " << (sizeof(B)/sizeof(*B)) << "\n";
    std::cout << "Size of C coefficient array: " << (sizeof(C)/sizeof(*C)) << "\n";
    // Print out initial conditions for the sim t=0
    std::cout << "Initial conditions for the nodes at time = 0 sec:\n";

    //***************************************//
    //******************SOLVE****************//
    //***************************************//

    // Itterate through the time steps
     while (timeCurr<=timeEnd)
    {
        outputFile<< "Time="<< timeCurr<<"," << "0";// far left node
        std::cout << "Values for the nodes at time = "<< timeCurr <<" sec:\n";
        std::cout << "0";
        for (int i=1;i<numPoints;i++)
        {
            std::cout << "," << F[i];
            outputFile << "," << F[i];
        }
        std::cout << " \n";
        outputFile << " \n";
        // Fill in the vlaues for the D array (as the b in the commonly known Ax=b equation)
        for (int i=0;i<numPoints-2;i++) // starts at node two F[1] 
        {
            D[i]=-F[i+1] -(alpha*deltat/(2*deltax*deltax)) * (F[i+2]-2*F[i+1]+F[i]);
        }

        C_star[0] = C[0] / B[0];
        D_star[0] = D[0] / B[0];

        // Do the forward sweep to reduce our matrix to a bidiagonal
        for (int i=1; i<numPoints-2;i++)
        {
            double m = 1.0/(B[i] - A[i] * C_star[i-1]);
            C_star[i] = C[i] * m;
            D_star[i] = (D[i] - A[i] * D_star[i-1]) * m;
        }

        // Finally we do a reverse sweep and set the new values of f
        for (int i=numPoints-2;i-- > 0;)
        {
            F[i+1] = D_star[i] - C_star[i] *D[i+1];
        }

    	timeCurr = timeCurr+deltat;
    }

    outputFile.close();
    return 0;
}
