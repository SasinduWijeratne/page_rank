#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include <math.h>
#include <string>
#include <cstring>
#include <limits>
#include <math.h>
#include "mpi.h"
using namespace std;

#define NUM_WORKERS 4
#define PREPROCESSING 0
/*mpic++ your_code_file.c
Execution

mpirun -np <no. of Processors> ./a.out
*/

const double DEFAULT_ALPHA = 0.85;
const double DEFAULT_CONVERGENCE = 0.00001;
const unsigned long DEFAULT_MAX_ITERATIONS = 100;
int trace = 1;

void read_partition_edges(char filename[], size_t* edgesDest){

    istream *infile;

    infile = new ifstream(filename);

    string line;

    int count = 0;
    while (getline(*infile, line)){
        int num = strtol(line.c_str(), NULL, 10);
        edgesDest[count] = num;
        count++;
    }
}

void read_partition_offset(char filename[], size_t* offsets){
    istream *infile;
    infile = new ifstream(filename);
    string line;
    int count = 0;
    while (getline(*infile, line)){
        int num = strtol(line.c_str(), NULL, 10);
        offsets[count] = num;
        count++;
    }
}

int main(int argc, char** argv) {
    int my_rank;  /* Identificador do processo */
    int proc_n;   /* Número de processos */
    int source;   /* Identificador do proc.origem */
    int dest;     /* Identificador do proc. destino */
    int tag = 50; /* Tag para as mensagens */
    double t1,t2; /* time stamps */
    int i, p, iteration;
    int MSIZE;
    FILE *f;

    MPI_Status status[500];  /* Status de retorno */
    MPI_Request request;
    MPI_Request requests[500];
    int recvd_count[500];
    int flag;

    MPI_Init (&argc , &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);
    p = proc_n - 1;

    if(my_rank == 0){

        string file_metadata = "partition0.metadata";

        size_t num_vertices;
        size_t tot_num_vertices;
        size_t num_edges;
        size_t s; // source

        FILE *fp = fopen(file_metadata.c_str(), "r");
        fscanf(fp, "%lu %lu %lu", &tot_num_vertices, &num_edges, &num_vertices);
        fclose(fp);

        printf("[Master] ");
        cout << tot_num_vertices << " " << num_edges << " " << num_vertices << endl;

        size_t each_num_vertices[NUM_WORKERS] = {0};

        for (int pid = 0; pid < NUM_WORKERS; pid++) {

            size_t left, right, num_vertices_per_worker = tot_num_vertices/NUM_WORKERS;
            left = num_vertices_per_worker*pid;
            right = (pid == NUM_WORKERS - 1) ? tot_num_vertices : num_vertices_per_worker*(pid+1);
            each_num_vertices[pid] = right - left;
        }

        size_t *edgesDest, *offsets, *recip_offsets;
        unsigned *labels;
        size_t *active_vertices;
        size_t length_active;

        labels = (unsigned *) calloc(tot_num_vertices + 1, sizeof(unsigned));
        active_vertices = (size_t *) calloc(tot_num_vertices , sizeof(size_t));
        
        cout << "Master alloc success\n";

        double convergence = DEFAULT_CONVERGENCE;
        unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;


        for (size_t k = 0; k < 1024; k++) {
            s = tot_num_vertices / 1024 * k;
            labels[s] = 1;
            active_vertices[k] = s;
        }
        length_active = 1024;

/*
        s = 0; // just an example
        labels[s] = 1;
        active_vertices[0] = s;
        length_active = 1;
*/

        int num_iterations = 0;
        int len = num_vertices;
        int num_rows = tot_num_vertices;

        while (num_iterations < max_iterations){
            double t1_iter, t2_iter;

            t1_iter = MPI_Wtime();

            if (length_active == 0) {
                for (i = 1; i < proc_n; i++)
                    MPI_Isend(active_vertices, 1, MPI_UNSIGNED_LONG, i, 51, MPI_COMM_WORLD, &requests[i]);
                for (i = 1; i < proc_n; i++)
                    MPI_Wait(&requests[i], &status[i]);   

                int count = 0;
                for (i = 0; i < tot_num_vertices; i++)
                    if (labels[i] == 0) {
                        count++;
                    }
                cout << "Unvisited: " << count;
                cout << endl;
                break;     
            }
/*
            cout << "length_active: " << length_active << endl;
            for (i = 0; i < length_active; i++)
                cout << active_vertices[i] << " ";
            cout << endl;
*/
            t1 = MPI_Wtime();
            for (i = 1; i < proc_n; i++)
                MPI_Isend(active_vertices, length_active, MPI_UNSIGNED_LONG, i, tag, MPI_COMM_WORLD, &requests[i]);
            for (i = 1; i < proc_n; i++) 
                MPI_Wait(&requests[i], &status[i]);
            
            t2 = MPI_Wtime();
            printf("[Master][Iter %d] Time for MPI send: %f\n", num_iterations, t2-t1);

            t1 = MPI_Wtime();
            for (i = 0; i < p; i++)
                MPI_Irecv(&active_vertices[i*each_num_vertices[0]], each_num_vertices[i], MPI_UNSIGNED_LONG, i+1, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[i]);
            for (i = 0; i < p; i++) {
                MPI_Wait(&requests[i], &status[i]);
                MPI_Get_count(&status[i], MPI_UNSIGNED_LONG, &recvd_count[i]);
            }
            t2 = MPI_Wtime();
            printf("[Master][Iter %d] Time for MPI recv: %f\n", num_iterations, t2-t1);

            t1 = MPI_Wtime();
            length_active = 0;
            for (i = 0; i < p; i++) {
                for (size_t j = 0; j < recvd_count[i]; j++) {
                    size_t vid = active_vertices[i*each_num_vertices[0]+j];
                    labels[vid] = num_iterations + 2; // when num_iterations == 0, we should label active vertices as 1
                    active_vertices[length_active] = vid;
                    length_active++;
                }
            }
            t2 = MPI_Wtime();
            printf("[Master][Iter %d] Time for Comp: %f\n", num_iterations, t2-t1);

            t2_iter = MPI_Wtime();
            printf("[Master][Iter %d] Iter time: %f\n", num_iterations, t2_iter-t1_iter);
            num_iterations++;            
        }
    }
    else {

        char file_metadata[30];
        char filename[20];
        char filename1[20];
        char filename2[20];

        sprintf(file_metadata, "partition%d.metadata", my_rank-1);
        sprintf(filename, "partition%d.txt", my_rank-1);
        sprintf(filename1, "offset%d.txt", my_rank-1);
        sprintf(filename2, "recipoffset.txt");

        size_t num_vertices;
        size_t tot_num_vertices;
        size_t num_edges;
        size_t left, right;

        FILE *fp = fopen(file_metadata, "r");
        fscanf(fp, "%lu %lu %lu", &tot_num_vertices, &num_edges, &num_vertices);
        fclose(fp);

        printf("[Worker %d] ", my_rank-1);
        cout << tot_num_vertices << " " << num_edges << " " << num_vertices << endl;

        left = (tot_num_vertices/NUM_WORKERS) * (my_rank - 1);
        right = (my_rank == NUM_WORKERS) ? tot_num_vertices : (tot_num_vertices/NUM_WORKERS)*my_rank;

        size_t *edgesDest, *offsets, *recip_offsets;
        size_t *active_vertices, *visited, *recvd_active_vertices;
        unsigned *labels;
        size_t length_active;

        edgesDest = (size_t *) malloc(num_edges * sizeof(size_t));
        offsets = (size_t *) malloc((tot_num_vertices) * sizeof(size_t));
        
        active_vertices = (size_t *) malloc(tot_num_vertices * sizeof(size_t));
        visited = (size_t *) calloc(num_vertices, sizeof(size_t));
        recvd_active_vertices = (size_t *) malloc(tot_num_vertices * sizeof(size_t));
        labels = (unsigned *) calloc(tot_num_vertices, sizeof(unsigned));

        printf("[Worker %d] alloc success\n", my_rank);

        read_partition_edges(filename,edgesDest);
        read_partition_offset(filename1,offsets);

        printf("[Worker %d] read success\n", my_rank);

        unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

        int num_iterations = 0;
        double diff_diff = 99;


        for (size_t k = 0; k < 1024; k++) {
            size_t s = tot_num_vertices / 1024 * k;
            labels[s] = 1;
            if (s >= left && s < right) {
                visited[s-left] = 1; 
            }
        }

/*
        size_t s = 0;
        labels[s] = 1;
        if (left == 0) {
            visited[0] = 1;
        }
*/
        while (num_iterations < max_iterations){
            double t1, t2;
            
            t1 = MPI_Wtime();
            MPI_Irecv(recvd_active_vertices, tot_num_vertices , MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[0]);
            MPI_Wait(&requests[0], &status[0]);
            
            MPI_Get_count(&status[0], MPI_UNSIGNED_LONG, &recvd_count[0]);
            t2 = MPI_Wtime();
            printf("[Node %d][Iter %d] Time for MPI recv: %f\n", my_rank, num_iterations, t2-t1);

            if (status[0].MPI_TAG == 51) {
                printf("[Node %d] Going to end\n", my_rank);
                break;
            }

            t1 = MPI_Wtime();

            for (i = 0; i < recvd_count[0]; i++)
                labels[recvd_active_vertices[i]] = num_iterations + 1;

            length_active = 0;
            for (i = 0; i < num_vertices; i++) {
                if (visited[i]) 
                    continue;
                for (size_t ci = offsets[i]; ci < offsets[i+1]; ci++) {
                    // check if edgesDest[ci] is in active list
                    if (labels[edgesDest[ci]] == num_iterations + 1) {
                        visited[i] = 1;
                        active_vertices[length_active] = i + left;
                        length_active++;
                        break;
                    }
                }
            }

            t2 = MPI_Wtime();
            printf("[Node %d][Iter %d] Time for Comp: %f\n", my_rank, num_iterations, t2-t1);

            t1 = MPI_Wtime();
            MPI_Isend(active_vertices, length_active, MPI_UNSIGNED_LONG, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status[0]);
            t2 = MPI_Wtime();
            printf("[Node %d][Iter %d] Time for MPI send: %f\n", my_rank, num_iterations, t2-t1);

            num_iterations++;

        }

    }

    MPI_Finalize();
    return 0;
}

