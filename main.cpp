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
#define SIZEV 640
#define MAXV 100
/*---------------------------------------------------------------------------*/
struct Edge
{
	int v1;
	int v2;
	int w;
	int s;
};
/*---------------------------------------------------------------------------*/
void readArr(int A[], int n, int rankn, int p, MPI_Comm comm);
void writeArr(int A[], int n, int rankn, int p, MPI_Comm comm);
void testPrintArr(int A[], int n, int rankn, int p, MPI_Comm comm);
void testMST(int A[], int n, int rankn, int p, MPI_Comm comm);
void printEdg(Edge* e, int n);

/*---------------------------------------------------------------------------*/

int main(int argc, char** argv)
{ 
	int rankn, p;
	int *A;
	int n = 14*4;
	MPI_Comm comm;
	double start, finish;

	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &p);
	MPI_Comm_rank(comm, &rankn);

	A = new int[n];
	readArr(A, n, rankn, p, comm);
	
	testPrintArr(A, n, rankn, p, comm);

	start = MPI_Wtime();
	testMST(A, n, rankn, p, comm);
	finish = MPI_Wtime();
	
	if (rankn == 0)
		cout << "Run time: " << finish-start << " seconds";
 
	/////////////////////////////////	
	delete[] A;
	MPI_Finalize();
	return 0;

}

void readArr(int A[], int n, int rankn, int p, MPI_Comm comm)
{
	int *temp = NULL;

	if (rankn == 0)
	{
		ifstream infile("graph.txt");
		int N,
			E;


		if(!infile.is_open())
			cout << "Cannot read file !!!" << endl;

		infile >> N;
		infile >> E;

		temp = new int[ E*4 ];


		for (int i = 0; i < E*4; i++)
		{
			infile >> temp[i++];
			infile >> temp[i++];
			infile >> temp[i++];
			temp[i] = 2;

		}
		infile.close();
	}
	
	MPI_Scatter(temp, n, MPI_INT, A, n, MPI_INT, 0, comm);
	
	if (rankn == 0)
	   delete[] temp;
}

void writeArr(int A[], int n, int rankn, int p, MPI_Comm comm)
{
	int* temp = NULL;
	
	if (rankn == 0)
	{
		ofstream outfile ("mst.txt");
		if (outfile.is_open())
		{
			temp = new int[n];

			MPI_Gather(A, n, MPI_INT, temp, n, MPI_INT, 0, comm);
			
			for (int i = 0; i < n; i++)
			{
				outfile << temp[i] << " ";
				outfile << temp[i++] << " ";
				outfile << temp[i++] << endl;
			}
			
			delete[] temp;
			outfile.close();
		}
	}
	else
	{
		MPI_Gather(A, n, MPI_INT, temp, n, MPI_INT, 0, comm);
	}

}

void testPrintArr(int A[], int n, int rankn, int p, MPI_Comm comm)
{
	int*       Temp;
	MPI_Status status;
	
	if (rankn == 0)
	{
		Temp = new int[n];
		cout << "N = " << n << endl;
		cout << "Rank: " <<  rankn << endl ;
		for (int i = 0; i < n; i++)
		{
			cout << "\t"<< A[i++] << " ";
			cout << A[i++] << " ";
			cout << A[i++] << " ";
			cout << A[i] << endl;
		}
		cout << endl;
	
		
		for (int j = 1; j < p; j++)
		{
			MPI_Recv(Temp, n, MPI_INT, j, 0, comm, &status);

			cout << "Rank: " <<  j << endl << "\t";
			for (int i = 0; i < n; i++)
				cout << Temp[i] << " ";
			cout << endl;
		}
		
		delete[] Temp;
	}
	else
	{
		MPI_Send(A, n, MPI_INT, 0, 0, comm);
	}
} 


void testMST(int A[], int n, int rankn, int p, MPI_Comm comm)
{
		
	int AM[SIZEV][MAXV];
	Edge* e = new Edge[n/4];
	Edge* mste = new Edge[n/4];

	for (int i=0; i < SIZEV; i++)
		for (int j=0; j < MAXV; j++)
			AM[i][j] = -1;

	if (rankn == 0)
	{
		int j = 0;
		for (int i = 0; i < n; i++)
		{
			e[j].v1 =  A[i++];
			e[j].v2 =  A[i++];
			e[j].w  =  A[i++];
			e[j].s  =  A[i];
			
			j++;
		}
		cout << "Rank: " <<  rankn << endl ;

		int iofe = 0;
		int c=0;
		for (c = 0; c < 9; c++)
		{
			int minv = 99999;
			cout << "Vertex : "<< c <<endl; 
			for (int i = 0; i < n/4; i++)
			{
				if ( (e[i].v1 == c ) || (e[i].v2 == c ) )
				{
					int x=0;
					for (x=0; AM[e[i].v1][x] >= 0; x++ );
					AM[ e[i].v1][x++] = e[i].v2; 
					AM[ e[i].v1][x] = e[i].w; 

					for (x=0; AM[e[i].v2][x] >= 0; x++ );
					AM[ e[i].v2][x++] = e[i].v2; 
					AM[ e[i].v2][x] = e[i].w; 


					if (e[i].w < minv) 
					{
						minv = e[i].w;
						iofe = i;
					}
//					cout << "\t Wight :" << e[i].w   << " - ";
//					cout << 		e[i].s   << " ";
//					cout << "\t" << e[i].v1 << " - " ;
//					cout << 		e[i].v2  << " ";
//					cout <<  endl;
				}

				
			}
//			cout << "\t The min wight is :" << e[iofe].w   << " - ";
//			cout << "\t" << e[iofe].v1 << " - " ;
//			cout << 		e[iofe].v2  << " ";
//			mste[c] = e[iofe];
//			cout <<  endl;
//			cout <<  endl;
		}
//
//	printEdg(mste, c);
//
		cout << endl;

	for (int i=0; i < SIZEV; i++)
	{
		if (AM[i][j] >=  0 )
		cout << " " << i << " -> " ;
		for (int j=0; j < MAXV; j++)
		{
			if (AM[i][j]>=0)cout << AM[i][j] << " ";
			if (AM[i][j] < 0 ) break;
		}
		if (AM[i][j] >=  0 )
		cout << endl;
	}
	cout << endl;

	}
	delete[] mste;
	delete[] e;
} 
void printEdg(Edge* e, int n)
{

	cout << "########################################" << endl;
	for (int i = 0; i < n; i++)
	{
	    cout << "\t Wight :" << e[i].w   << " - ";
	    cout << 		e[i].s   << " ";
	    cout << "\t" << e[i].v1 << " - " ;
	    cout << 		e[i].v2  << " ";
	    cout <<  endl;
	}
	cout <<  endl;
}
