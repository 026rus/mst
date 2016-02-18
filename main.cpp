/* 
 * File:   main.cpp
 * Author: Vitaly
 *
 * Created on Tue Feb  9 10:00:39 2016
 *
 */

#include <iostream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <mpi.h>

using namespace std;

/*
 * 
 */
/*---------------------------------------------------------------------------*/
bool debug = false;
/*---------------------------------------------------------------------------*/
int findp(int *subsets, int i);
void Union(int * subsets, int x, int y, int &);
void MST(int ** edge, int V, int E, int * &, int &,  int rankn, int p, MPI_Comm comm);
void readArr(int** &A, int &V, int &E, int rankn, int p, MPI_Comm comm);
void writeArr(int A[], int n);
void flatit(int **&A, int * &Fl, int E);
void unflatit(int *&A, int ** &,  int E);
void ftf(int ** &A, int *&, int &, int &, int rankn, int p, MPI_Comm comm);
void printA(int *A, int N);
/*---------------------------------------------------------------------------*/

int main(int argc, char** argv)
{ 
	int rankn=0,
		p=0;
	int **A;
	int * MSTarr;
	int MSTcount;
	int *Fl;
	int V=0,
		E=0;
	MPI_Comm comm;
	double start, finish;

	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &p);
	MPI_Comm_rank(comm, &rankn);

	if (!rankn)
		readArr(A, V, E, rankn, p, comm);

	start = MPI_Wtime();

	ftf(A, Fl, E, V ,rankn, p, comm);
	
 	MST(A, V, E, MSTarr, MSTcount, rankn, p, comm);

	finish = MPI_Wtime();
	if (rankn == 0)
		cout << "Run time: " << finish-start << " seconds";
	if (!rankn)
		writeArr(MSTarr, MSTcount);
	/////////////////////////////////	
/*
	if (A != NULL)
	{
		for (int i=0; i<E; ++i)
			delete[] A[i];
		delete[] A;

	}
	*/
	MPI_Finalize();
	return 0;

}
/*************************************************************************************************/
/* Sending array of edges to all procesors */
void ftf(int **& A, int *& Fl, int &E, int &V, int rankn, int p, MPI_Comm comm)
{
	MPI_Bcast(&E, 1, MPI_INT, 0, comm);
	MPI_Bcast(&V, 1, MPI_INT, 0, comm);

	if(rankn)
		Fl = new int[E*3];

	if(rankn)
	{
		A = new int*[E];
		for (int i=0; i< E; ++i)
			A[i] = new int [3];
	}
	
	if (!rankn)
		flatit(A, Fl, E);

	MPI_Bcast(Fl, E*3, MPI_INT, 0, comm);

	if(rankn)
		unflatit(Fl, A, E);

	MPI_Barrier(comm);
}
/**************************************************************************************************/
/* converting 2D arr in to 1D arr 
 * to make brotcast esier
 */
void flatit(int **&A, int * & temp, int E)
{
	temp = new int[E*3];

	int c=0;
	for (int i=0; i< E*3; i+=3)
	{
		temp[i] 	=	A[c][0];
		temp[i+1] 	=	A[c][1];
		temp[i+2] 	=	A[c++][2];
	}
}
/**************************************************************************************************/
/* Converting 1D array int to 2D  */
void unflatit(int *&T, int ** &A, int E)
{
	int c=0;
	for (int i=0; i< E; ++i)
	{
		A[i][0] = 	T[c++];
		A[i][1] =	T[c++];
		A[i][2] =	T[c++];
	}
}
/**************************************************************************************************/
void writeArr(int A[], int n)
{
	ofstream outfile ("mst.txt");
	if (outfile.is_open())
	{
		for (int i = 0; i < n; i++)
		{
			outfile << A[i++] << " - ";
			outfile << A[i] << endl;
		}
		
		outfile.close();
	}

}
/**************************************************************************************************/
/* reading array from the file */
void readArr(int** &A, int &V, int &E, int rankn, int p, MPI_Comm comm)
{

	ifstream infile("graph.txt");


	if(!infile.is_open())
		cout << "Cannot read file !!!" << endl;

	infile >> V;
	infile >> E;

	A = new int*[ E ];
	for (int i=0; i < E;  ++i)
		A[i] = new int[3];


	for (int c = 0; c < E; ++c)
	{
	
		infile >> A[c][0];	// v1
		infile >> A[c][1];	// v2
		infile >> A[c][2];	// wight
	
	}
	infile.close();
}
/**************************************************************************************************/
/* The MST function :) */
void MST(int ** A, int V, int E, int * &TheMST, int &MSTcount, int rankn, int p, MPI_Comm comm)
{
	int sizep	= V/p;
	int startp 	= sizep*rankn;
	int endp 	= startp + sizep;
	TheMST = new int[V*2];
	int * minset = new int[V];
	int * minset2 = new int[sizep];
	int * pointlist = new int[V];

    for (int v = 0; v < V; ++v)
		pointlist[v] = v;
 
    int numTrees = V;
    int MSTweight = 0;	// For debug
    while (numTrees > 1)
    {
		// set minset to empty
		for (int v = 0; v < V; ++v)
				minset[v] = -1;

		for (int i=0; i<sizep; ++i)
			minset2[i] = -1;
		/*----------------------------------------------*/
        for (int i=0; i<E; i++)
        {
            int set1 = findp(pointlist, A[i][0]);
            int set2 = findp(pointlist, A[i][1]);
 
            if (set1 != set2)
			{
				if (set1 >= startp && set1 < endp)
				{
					if (minset2[set1%sizep] == -1 || A[minset2[set1%sizep]][2] > A[i][2])
						minset2[set1%sizep] = i;
				}	
				if (set2 >= startp && set2 < endp)
				{
					if (minset2[set2%sizep] == -1 || A[minset2[set2%sizep]][2] > A[i][2])
						minset2[set2%sizep] = i;
				}
			}
        }
		/*----------------------------------------------*/
		MPI_Gather(minset2, sizep, MPI_INT, minset, sizep, MPI_INT, 0, comm);
		MPI_Barrier(comm);
		/*----------------------------------------------*/
 		if(!rankn)
		{

        	for (int i=0; i<V; i++)
        	{
        	    if (minset[i] != -1)
        	    {
        	        int set1 = findp(pointlist, A[minset[i]][0]);
        	        int set2 = findp(pointlist, A[minset[i]][1]);
 
        	        if (set1 != set2)
					{
						if (debug) MSTweight += A[minset[i]][2];
						TheMST[MSTcount++] = A[minset[i]][0]; 
						TheMST[MSTcount++] = A[minset[i]][1]; 
						Union(pointlist, set1, set2, numTrees);
					}
        	    }
        	}
		}
		MPI_Bcast(&numTrees, 1, MPI_INT, 0, comm);
		MPI_Bcast(pointlist, V, MPI_INT, 0, comm);
    }
 
    if (debug) cout << "Weight of MST is " << MSTweight << endl;
	if(!rankn)

	if (debug)	printA(TheMST, MSTcount);
	delete [] pointlist;
}
/**************************************************************************************************/
/* Find the perent of node */
int findp(int * s, int i)
{
    if (s[i] != i)
		s[i] = findp(s, s[i]);
 
    return s[i];
}
/**************************************************************************************************/
/* Union 2 vertises */
void Union(int* s, int x, int y, int & numTrees)
{
    int xroot = findp(s, x);
    int yroot = findp(s, y);
 
    s[yroot] = xroot;

	numTrees--;
}
/**************************************************************************************************/
/* print array for testing only  */
void printA(int *e, int E)
{
	if (e == NULL) return;

	for (int i = 0; i < E; i++)
	{
	    cout << 		e[i++]  << " - ";
	    cout << 		e[i] ;
	    cout <<  endl;
	}
	cout <<  endl;
}
