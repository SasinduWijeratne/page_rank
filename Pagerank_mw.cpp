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

#define NUM_WORKERS 16
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

void read_partition_write(size_t *offsets, size_t *edgesDest, int tot_num_vertices){
    FILE *f;
    char filename[20];
    for (int pid = 0; pid < NUM_WORKERS; pid++) {
        sprintf(filename, "partition%d.txt", pid);
        f = fopen(filename, "w");
        size_t left, right, num_vertices_per_worker = tot_num_vertices/NUM_WORKERS;
        left = num_vertices_per_worker*pid;
        right = (pid == NUM_WORKERS - 1) ? tot_num_vertices : num_vertices_per_worker*(pid+1);

        for (size_t k = offsets[left]; k < offsets[right]; k++)
            fprintf(f, "%lu\n", edgesDest[k]);
        fclose(f);
    }
}

void partition_meta(size_t tot_vertices, size_t *offsets){
    FILE *f;
    char filename[20];

    for (int pid = 0; pid < NUM_WORKERS; pid++) {
        sprintf(filename, "partition%d.metadata", pid);
        f = fopen(filename, "w");
        size_t left, right, num_vertices_per_worker = tot_vertices/NUM_WORKERS;
        left = num_vertices_per_worker*pid;
        right = (pid == NUM_WORKERS - 1) ? tot_vertices : num_vertices_per_worker*(pid+1);
        fprintf(f, "%lu %lu %lu\n", tot_vertices, offsets[right] - offsets[left], right - left);
        fclose(f);
    }
}

void read_offset_write(size_t *offsets, int tot_size){ 
    FILE *f;
    char filename[20];

    size_t num_vertices_per_worker = tot_size / NUM_WORKERS;

    for (int pid = 0; pid < NUM_WORKERS; pid++) {
        sprintf(filename, "offset%d.txt", pid);
        f = fopen(filename, "w");
        size_t left, right;
        left = num_vertices_per_worker*pid;
        right = (pid == NUM_WORKERS - 1) ? tot_size : num_vertices_per_worker*(pid+1);
        for (int i = left; i <= right; i++)
            fprintf(f, "%zu\n", offsets[i] - offsets[left]);  
        fclose(f);
    }
}

void read_recipoffset_write( size_t *offsets, int tot_size){ 
    FILE *f;
    char filename[20];
    sprintf(filename, "recipoffset.txt");
    f = fopen(filename, "w");
    for(int i = 0; i <= tot_size; i++ ){
        fprintf(f, "%zu\n", offsets[i]);
    }
    fclose(f);
}

void read_matrix(const string &filename, size_t *edgesDest, size_t *offsets, size_t *recip_offsets, int num_vertices, int num_edges){
    istream *infile;
    infile = new ifstream(filename.c_str());
    string line;

    size_t* temp_offsets;
    temp_offsets = (size_t *) calloc((num_vertices) , sizeof(size_t));

    size_t count = 0;
    size_t prev = -1;
    while (getline(*infile, line)){
        stringstream lineStream(line);
        size_t src, dest;
        lineStream >> src >> dest;
        offsets[dest+1]++;
        if (prev != src) {
            prev++;
            for (; prev < src; prev++)
                recip_offsets[prev] = count;
            recip_offsets[src] = count;
            prev = src;
        }
        count++;
    }
    prev++;
    for (; prev < num_vertices; prev++)
        recip_offsets[prev] = count;

    offsets[0] = 0;
    for(int kk = 1; kk < num_vertices; kk++)
        offsets[kk] = offsets[kk-1] + offsets[kk];

    istream *infile2;
    infile2 = new ifstream(filename.c_str());
    string line2;

    while (getline(*infile2, line2)){
        stringstream lineStream(line2);
        size_t src, dest;
        lineStream >> src >> dest;
        edgesDest[offsets[dest] + temp_offsets[dest]] = src;
        temp_offsets[dest]++;
    }
    free(temp_offsets);
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

        size_t each_num_vertices[NUM_WORKERS] = {0};

        for (int pid = 0; pid < NUM_WORKERS; pid++) {

            size_t left, right, num_vertices_per_worker = tot_num_vertices/NUM_WORKERS;
            left = num_vertices_per_worker*pid;
            right = (pid == NUM_WORKERS - 1) ? tot_num_vertices : num_vertices_per_worker*(pid+1);
            each_num_vertices[pid] = right - left;
        }

        double *pr, *old_pr;
        size_t *edgesDest, *offsets, *recip_offsets;

        pr = (double *) calloc(tot_num_vertices +1, sizeof(double));
        old_pr = (double *) calloc(tot_num_vertices , sizeof(double));

        double convergence = DEFAULT_CONVERGENCE;
        unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

        double sum_pr;
        double dangling_pr;
        double alpha = DEFAULT_ALPHA;

        sum_pr = 0;
        dangling_pr = 0;
        pr[0] = 1;

        int num_iterations = 0;
        double diff = -1;
        double diff_prev = -1;
        int len = num_vertices;
        int num_rows = tot_num_vertices;
        double diff_diff = 99;

        while ((diff_diff > convergence) && num_iterations < max_iterations){
            double t1_iter, t2_iter;
            t1_iter = MPI_Wtime();

            diff_prev = diff;
            diff = 0;
            for(i = 0; i < tot_num_vertices; i++)
                diff += fabs(pr[i] - old_pr[i]);
            diff_diff = fabs(diff - diff_prev);

            printf("[Iter %d] Diff: %f\n", num_iterations, diff);
            
            // printf("[ID 0] diff: %f\n",diff_diff);
            pr[tot_num_vertices] = diff_diff;

            t1 = MPI_Wtime();
            for (i = 1; i < proc_n; i++)
                MPI_Isend(&pr[0], tot_num_vertices + 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &requests[i]);
            for (i = 1; i < proc_n; i++)
                MPI_Wait(&requests[i], &status);
            t2 = MPI_Wtime();
            printf("[Master] Time for MPI send: %f\n", t2-t1);
            for(int ii = 0; ii < tot_num_vertices; ii++)
                old_pr[ii] = pr[ii];

            t1 = MPI_Wtime();
            for (i = 0; i < p; i++)
                MPI_Irecv(&pr[i*each_num_vertices[0]], each_num_vertices[i], MPI_DOUBLE, i+1, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[i]);
            for (i = 0; i < p; i++)
                MPI_Wait(&requests[i], &status);
            t2 = MPI_Wtime();
            printf("[Master] Time for MPI recv: %f\n", t2-t1);

            num_iterations++;
            
            t2_iter = MPI_Wtime();
            printf("[Master] Iter time: %f\n", t2_iter-t1_iter);

        }

    }
    else{

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

        FILE *fp = fopen(file_metadata, "r");
        fscanf(fp, "%lu %lu %lu", &tot_num_vertices, &num_edges, &num_vertices);
        fclose(fp);

        double *pr, *new_pr;
        size_t *edgesDest, *offsets, *recip_offsets;

        edgesDest = (size_t *) malloc(num_edges * sizeof(size_t));
        offsets = (size_t *) malloc((tot_num_vertices) * sizeof(size_t));
        recip_offsets = (size_t *) malloc((tot_num_vertices) * sizeof(size_t));

        read_partition_edges(filename,edgesDest);
        read_partition_offset(filename1,offsets);
        read_partition_offset(filename2,recip_offsets);

        pr = (double *) calloc(tot_num_vertices+1 , sizeof(double));
        new_pr = (double *) calloc(num_vertices , sizeof(double));

        double convergence = DEFAULT_CONVERGENCE;
        unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;


        double dangling_pr;
        double alpha = DEFAULT_ALPHA;

        dangling_pr = 0;
        pr[0] = 1;

        int num_iterations = 0;
        double diff_diff = 99;

        while ((diff_diff > convergence) && num_iterations < max_iterations){

            double t1_iter, t2_iter;
            t1_iter = MPI_Wtime();
            
            double t1, t2;
            t1 = MPI_Wtime();
            MPI_Irecv(&pr[0], tot_num_vertices+1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            t2 = MPI_Wtime();
            printf("[Node %d][Iter %d] Time for MPI recv: %f\n", my_rank, num_iterations, t2-t1);

            t1 = MPI_Wtime();

            diff_diff = pr[tot_num_vertices];
            // printf("[ID %d] diff: %f\n",my_rank,diff_diff);

            for (i = 0; i < num_vertices; i++) {
                /* The corresponding element of the H multiplication */
                double h = 0.0;
                for (size_t ci = offsets[i]; ci < offsets[i+1]; ci++) {
                    /* The current element of the H vector */
                    double h_v = (recip_offsets[edgesDest[ci]+1]-recip_offsets[edgesDest[ci]])
                        ? 1.0 / (recip_offsets[edgesDest[ci]+1]-recip_offsets[edgesDest[ci]])
                        : 0.0;
                    h += h_v * pr[edgesDest[ci]];
                }
                h *= alpha;
                new_pr[i] = h + (1-alpha)/tot_num_vertices;
            }
            
            t2 = MPI_Wtime();
            printf("[Node %d][Iter %d] Time for Comp: %f\n", my_rank, num_iterations, t2-t1);

            t1 = MPI_Wtime();
            MPI_Isend(&new_pr[0], num_vertices, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            t2 = MPI_Wtime();
            printf("[Node %d][Iter %d] Time for MPI send: %f\n", my_rank, num_iterations, t2-t1);

            num_iterations++;
            
            t2_iter = MPI_Wtime();
            printf("[Node %d][Iter %d] Iter time: %f\n", my_rank, num_iterations, t2_iter-t1_iter);

        }

    }

    MPI_Finalize();
    return 0;
}

