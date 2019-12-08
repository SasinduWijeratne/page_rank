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

    MPI_Status status;  /* Status de retorno */
    MPI_Request request;
    MPI_Request requests[500];
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

        FILE *fp = fopen(file_metadata.c_str(), "r");
        fscanf(fp, "%lu %lu %lu", &tot_num_vertices, &num_edges, &num_vertices);
        fclose(fp);

        cout << tot_num_vertices << " " << num_edges << " " << num_vertices << endl;

        size_t each_num_vertices[NUM_WORKERS] = {0};

        for (int pid = 0; pid < NUM_WORKERS; pid++) {
            size_t left, right, num_vertices_per_worker = tot_num_vertices/NUM_WORKERS;
            left = num_vertices_per_worker*pid;
            right = (pid == NUM_WORKERS - 1) ? tot_num_vertices : num_vertices_per_worker*(pid+1);
            each_num_vertices[pid] = right - left;
        }

        size_t *pr, *old_pr;
        size_t *edgesDest, *offsets, *recip_offsets;

        pr = (size_t *) malloc((tot_num_vertices +1)* sizeof(size_t));
        old_pr = (size_t *) malloc((tot_num_vertices +1)* sizeof(size_t));

        unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

        for (size_t k = 0; k < tot_num_vertices; k++)
            pr[k] = tot_num_vertices;

        int num_iterations = 0;

        int len = num_vertices;
        int num_rows = tot_num_vertices;
        int hasUpdate = 1;


        cout << "pr[0]: " << pr[0] << endl;
        cout << pr[46372] << " " << pr[84294] << endl;

        while (hasUpdate && num_iterations < max_iterations){
            double t1_iter, t2_iter;
            t1_iter = MPI_Wtime();
            hasUpdate = 0;
            for(int ii = 0; ii < tot_num_vertices; ii++) {
                if (old_pr[ii] != pr[ii]) {
                    hasUpdate = 1;
                    break;
                }
            }
            pr[tot_num_vertices] = hasUpdate;
            t1 = MPI_Wtime();
            for (i = 1; i < proc_n; i++)
                MPI_Isend(&pr[0], tot_num_vertices + 1, MPI_UNSIGNED_LONG, i, tag, MPI_COMM_WORLD, &requests[i]);
            for (i = 1; i < proc_n; i++)
                MPI_Wait(&requests[i], &status);
            t2 = MPI_Wtime();
            printf("[Master][Round 1] Time for MPI send: %f\n", t2-t1);
            for(int ii = 0; ii < tot_num_vertices; ii++)
                old_pr[ii] = pr[ii];

            t1 = MPI_Wtime();
            for (i = 0; i < p; i++)
                MPI_Irecv(&pr[i*each_num_vertices[0]], each_num_vertices[i], MPI_UNSIGNED_LONG, i+1, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[i]);
            for (i = 0; i < p; i++)
                MPI_Wait(&requests[i], &status);
            t2 = MPI_Wtime();
            printf("[Master][Round 1] Time for MPI recv: %f\n", t2-t1);

            num_iterations++;

            int count = 0;
            for (int kk = 0; kk < tot_num_vertices; kk++)
                if (pr[kk] != 0) {
                    count++;
                } else {
                    if (kk-count <= 10)
                        cout << kk << " ";
                }
            cout << endl;
            cout << "Non-zero: " << count <<endl;
            
            t2_iter = MPI_Wtime();
            printf("[Master][Round 1] Iter time: %f\n", t2_iter-t1_iter);
        }

        num_iterations = 0;
        hasUpdate = 1;
        while ((hasUpdate || num_iterations < 3) && num_iterations < max_iterations){
            double t1_iter, t2_iter;
            t1_iter = MPI_Wtime();
            hasUpdate = 0;
            for(int ii = 0; ii < tot_num_vertices; ii++) {
                if (old_pr[ii] != pr[ii]) {
                    hasUpdate = 1;
                    break;
                }
            }
            pr[tot_num_vertices] = hasUpdate;
            t1 = MPI_Wtime();
            for (i = 1; i < proc_n; i++)
                MPI_Isend(&pr[0], tot_num_vertices + 1, MPI_UNSIGNED_LONG, i, tag, MPI_COMM_WORLD, &requests[i]);
            for (i = 1; i < proc_n; i++)
                MPI_Wait(&requests[i], &status);
            t2 = MPI_Wtime();
            printf("[Master][Round 2] Time for MPI send: %f\n", t2-t1);
            for(int ii = 0; ii < tot_num_vertices; ii++)
                old_pr[ii] = pr[ii];

            t1 = MPI_Wtime();
            for (i = 0; i < p; i++)
                MPI_Irecv(&pr[i*each_num_vertices[0]], each_num_vertices[i], MPI_UNSIGNED_LONG, i+1, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[i]);
            for (i = 0; i < p; i++)
                MPI_Wait(&requests[i], &status);
            t2 = MPI_Wtime();
            printf("[Master][Round 2] Time for MPI recv: %f\n", t2-t1);

            int count = 0;
            for (int kk = 0; kk < tot_num_vertices; kk++)
                if (pr[kk] != 0) {
                    count++;
                } else {
                    if (kk-count <= 10)
                        cout << kk << " ";
                }
            cout << endl;
            cout << "Non-zero: " << count <<endl;
            
            num_iterations++;
            
            t2_iter = MPI_Wtime();
            printf("[Master][Round 2] Iter time: %f\n", t2_iter-t1_iter);
        }


    }
    else{

        char file_metadata[30];
        char file_metadata2[30];
        char filename[30];
        char filename1[30];
        char filename2[30];
        char filename3[30];

        sprintf(file_metadata, "partition%d.metadata", my_rank-1);
        sprintf(file_metadata2, "gpoppartition%d.metadata", my_rank-1);
        sprintf(filename, "partition%d.txt", my_rank-1);
        sprintf(filename1, "offset%d.txt", my_rank-1);
        sprintf(filename2, "gpoppartition%d.txt", my_rank-1);
        sprintf(filename3, "gpopoffset%d.txt", my_rank-1);

        size_t num_vertices;
        size_t tot_num_vertices;
        size_t num_edges, num_edges2;

        FILE *fp = fopen(file_metadata, "r");
        fscanf(fp, "%lu %lu %lu", &tot_num_vertices, &num_edges2, &num_vertices);
        fclose(fp);

        cout << "[Worker " << my_rank << "] " << tot_num_vertices << " " << num_edges2 << " " << num_vertices << endl;

        fp = fopen(file_metadata2, "r");
        fscanf(fp, "%lu %lu %lu", &tot_num_vertices, &num_edges, &num_vertices);
        fclose(fp);

        cout << "[Worker " << my_rank << "] " <<  tot_num_vertices << " " << num_edges << " " << num_vertices << endl;

        size_t *pr, *new_pr, *rank;
        size_t *edgesDest, *offsets, *edgesDest2, *offsets2;

        edgesDest = (size_t *) malloc(num_edges * sizeof(size_t));
        offsets = (size_t *) malloc((tot_num_vertices) * sizeof(size_t));
        edgesDest2 = (size_t *) malloc(num_edges2 * sizeof(size_t));
        offsets2 = (size_t *) malloc((tot_num_vertices) * sizeof(size_t));

        read_partition_edges(filename,edgesDest2);
        read_partition_offset(filename1,offsets2);
        read_partition_edges(filename2,edgesDest);
        read_partition_offset(filename3,offsets);

        pr = (size_t *) malloc((tot_num_vertices+1)*sizeof(size_t));
        new_pr = (size_t *) malloc(num_vertices*sizeof(size_t));
        rank = (size_t *) malloc(num_vertices*sizeof(size_t));

        size_t left, num_vertices_per_worker = tot_num_vertices/NUM_WORKERS;
        left = num_vertices_per_worker*(my_rank-1);
        
        for (size_t v = 0; v < num_vertices; v++)
            rank[v] = v + left;

        unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

        int num_iterations = 0;
        int hasUpdate = 1;

        while (hasUpdate && num_iterations < max_iterations){
            double t1_iter, t2_iter;
            t1_iter = MPI_Wtime();
            
            double t1, t2;
            t1 = MPI_Wtime();
            MPI_Irecv(&pr[0], tot_num_vertices+1, MPI_UNSIGNED_LONG, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            t2 = MPI_Wtime();
            printf("[Node %d][Round 1][Iter %d] Time for MPI recv: %f\n", my_rank, num_iterations, t2-t1);

            t1 = MPI_Wtime();

            hasUpdate = pr[tot_num_vertices];
            for (size_t v = 0; v < num_vertices; v++) {
                for (size_t ci = offsets2[v]; ci < offsets2[v+1]; ci++) {
                    size_t u = edgesDest2[ci];
                    if (pr[u] < rank[v]) 
                        rank[v] = pr[u];
                }
                new_pr[v] = rank[v];
            }
            t2 = MPI_Wtime();
            printf("[Node %d][Round 1][Iter %d] Time for Comp: %f\n", my_rank, num_iterations, t2-t1);

            t1 = MPI_Wtime();
            MPI_Isend(&new_pr[0], num_vertices, MPI_UNSIGNED_LONG, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            t2 = MPI_Wtime();
            printf("[Node %d][Round 1][Iter %d] Time for MPI send: %f\n", my_rank, num_iterations, t2-t1);

            num_iterations++;
            
            t2_iter = MPI_Wtime();
            printf("[Node %d][Round 1][Iter %d] Iter time: %f\n", my_rank, num_iterations, t2_iter-t1_iter);
        }

        num_iterations = 0;
        hasUpdate = 1;
        while ((hasUpdate || num_iterations < 3) && num_iterations < max_iterations){
            double t1_iter, t2_iter;
            t1_iter = MPI_Wtime();
            
            double t1, t2;
            t1 = MPI_Wtime();
            MPI_Irecv(&pr[0], tot_num_vertices+1, MPI_UNSIGNED_LONG, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            t2 = MPI_Wtime();
            printf("[Node %d][Round 2][Iter %d] Time for MPI recv: %f\n", my_rank, num_iterations, t2-t1);

            t1 = MPI_Wtime();
            hasUpdate = pr[tot_num_vertices];
            for (size_t v = 0; v < num_vertices; v++) {
                for (size_t ci = offsets[v]; ci < offsets[v+1]; ci++) {
                    size_t u = edgesDest[ci];
                    if (pr[u] < rank[v]) 
                        rank[v] = pr[u];
                }
                new_pr[v] = rank[v];
            }
            t2 = MPI_Wtime();
            printf("[Node %d][Round 2][Iter %d] Time for Comp: %f\n", my_rank, num_iterations, t2-t1);

            t1 = MPI_Wtime();
            MPI_Isend(&new_pr[0], num_vertices, MPI_UNSIGNED_LONG, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            t2 = MPI_Wtime();
            printf("[Node %d][Round 2][Iter %d] Time for MPI send: %f\n", my_rank, num_iterations, t2-t1);

            num_iterations++;
            
            t2_iter = MPI_Wtime();
            printf("[Node %d][Round 2][Iter %d] Iter time: %f\n", my_rank, num_iterations, t2_iter-t1_iter);

        }


    }

    MPI_Finalize();
    return 0;
}

