#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <chrono>

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
	cout << "# ";
	for (int i = 0; i < argc; i++) cout << argv[i] << " ";
	cout << endl;

    MPI_Init(&argc, &argv);
 
    int my_rank; int np;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // get rank of my process
	MPI_Comm_size(MPI_COMM_WORLD, &np);	// get total number of processes
 
	int* buffer = new int[np]();  // buffer to hold all relations from a loop iteration
	srand(time(NULL)^my_rank);  // seed RNG from time, rank
	int n = 5 + (rand() % 5);  // choose a random loop length in {5,6,7,8,9}

	MPI_Comm comm1 = MPI_COMM_WORLD;
	MPI_Comm comm2;

	for (int i = 0; i < n; i++) {
		MPI_Comm_rank(comm1, &my_rank);
		MPI_Comm_size(comm1, &np);

		int rel = rand() % 30000;  // we have a random integer rel representing a relation
		// now send this rank's rel to all ranks
		MPI_Allgather(&rel, 1, MPI_INT, buffer, 1, MPI_INT, comm1);
	
		cout << n << "," << my_rank << ":";
		for (int j = 0; j < np; j++) {
			cout << buffer[j];
			if (j < np-1) cout << ",";
		}
		cout << endl;

		// update communicator if at end of ragged loop
		int col = 0; int key = 1;
		if (i == n-1) col = 0;
		else col = 1;
		
		MPI_Comm_split(comm1, col, key, &comm2);
		comm1 = comm2;
	}
 
    MPI_Finalize();
 
    return EXIT_SUCCESS;
}
