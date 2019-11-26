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
#define NUM_WORKERS 128
#define PREPROCESSING 0
/*mpic++ your_code_file.c
Execution

mpirun -np <no. of Processors> ./a.out
*/

const double DEFAULT_ALPHA = 0.85;
const double DEFAULT_CONVERGENCE = 0.001;
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

#if PREPROCESSING == 0
void ring_pagerank(int id, MPI_Status status, int proc_n, int tag){

    double convergence = DEFAULT_CONVERGENCE;
    unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

    double sum_pr;
    double dangling_pr;
    double alpha = DEFAULT_ALPHA;
    int i;

    sum_pr = 0;
    dangling_pr = 0;

    size_t num_vertices;
    size_t tot_num_vertices;
    size_t num_edges;

    MPI_Request request;

    char file_metadata[30];
    char filename[20];
    char filename1[20];
    char filename2[20];
    sprintf(file_metadata, "partition%d.metadata", id);
    sprintf(filename, "partition%d.txt", id);
    sprintf(filename1, "offset%d.txt", id);
    sprintf(filename2, "recipoffset.txt");

    FILE *fp = fopen(file_metadata, "r");
    fscanf(fp, "%lu %lu %lu", &tot_num_vertices, &num_edges, &num_vertices);
    fclose(fp);

    double  *pr, *old_pr;
    size_t *edgesDest, *offsets, *recipoffsets;
    size_t left, right;

    edgesDest = (size_t *) calloc(num_edges , sizeof(size_t));
    offsets = (size_t *) calloc((tot_num_vertices) , sizeof(size_t));
    recipoffsets = (size_t *) calloc((tot_num_vertices) , sizeof(size_t));

    pr = (double *) calloc(tot_num_vertices/NUM_WORKERS + NUM_WORKERS, sizeof(double));
    old_pr = (double *) calloc(tot_num_vertices, sizeof(double));

    double *pr1;
    pr1 = (double *) malloc((tot_num_vertices/NUM_WORKERS + NUM_WORKERS) * sizeof(double));

    read_partition_edges(filename,edgesDest);
    read_partition_offset(filename1,offsets);
    read_partition_offset(filename2,recipoffsets);


    if (id == 0)
        pr[0] = 1;

    int num_iterations = 0;
    double diff = 1;
    double diff_prev = 999;
    int exit_flag = 0;
    // (fabs(diff - diff_prev) > convergence)
    struct timespec start_iter, stop_iter, start, stop;
    double time_iter, time;

    while ((num_iterations < max_iterations) && diff > DEFAULT_CONVERGENCE && fabs(diff - diff_prev) > 0.00001 ) {
        printf("[ID %d] Iter: %d Diff: %f \n", id, num_iterations, diff);

        if( clock_gettime( CLOCK_REALTIME, &start_iter) == -1 ) { perror( "clock gettime" );}
        
        if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
        left = (tot_num_vertices/NUM_WORKERS)*id;

        if (num_iterations >= 1) {
            if (num_iterations >= 2)
                diff_prev = diff;
            diff = 0;
            for (size_t k = 0; k < num_vertices; k++)
                diff += fabs(old_pr[k + left] - pr[k]);
        }

        if (num_iterations == 0)
            old_pr[0] = 1;

        for (i = 0; i < num_vertices; i++)
            pr[i] = (1-alpha)/tot_num_vertices;

        for (i = 0; i < num_vertices; i++) {
            double h = 0.0;
            for (size_t ci = offsets[i]; ci < offsets[i+1]; ci++) {
                double h_v = (recipoffsets[edgesDest[ci]+1]-recipoffsets[edgesDest[ci]])
                                ? 1.0 / (recipoffsets[edgesDest[ci]+1]-recipoffsets[edgesDest[ci]])
                                : 0.0;
                h += h_v * old_pr[edgesDest[ci]];
            }
            h *= alpha;
            pr[i] += h;
        }
        pr[tot_num_vertices/NUM_WORKERS + NUM_WORKERS - 1] = diff;

        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}       
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("[ID %d][Iter %d] Comp Time: %f sec\n", id, num_iterations, time);

        diff = 0;
        
        if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
        for(int sas_i = 0; sas_i < proc_n; sas_i++){

            int ring_partition_id = (id + proc_n - sas_i)%proc_n;
            left = (tot_num_vertices/NUM_WORKERS)*ring_partition_id;
            right = (ring_partition_id == NUM_WORKERS - 1) ? tot_num_vertices : (tot_num_vertices/NUM_WORKERS)*(ring_partition_id+1);
            int current_partition_size = right - left;

            // if sas_i != proc_n
            if (id%2 == 0) {
                for (size_t k = 0; k < current_partition_size; k++)
                    old_pr[left + k] = pr[k];
                diff += pr[tot_num_vertices/NUM_WORKERS + NUM_WORKERS - 1];

                MPI_Isend(&pr[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_DOUBLE, (id+1)%proc_n, tag, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                int use_id = (id == 0)? (proc_n-1): (id-1)%proc_n;

                MPI_Irecv(&pr[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_DOUBLE, use_id, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
            }
            else {

                MPI_Irecv(&pr1[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_DOUBLE,  (id-1)%proc_n, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);

                for (size_t k = 0; k < current_partition_size; k++)
                    old_pr[left + k] = pr[k];
                diff += pr[tot_num_vertices/NUM_WORKERS + NUM_WORKERS - 1];

                MPI_Isend(&pr[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_DOUBLE, (id+1)%proc_n, tag, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);

                double * swap = pr;
                pr = pr1;
                pr1 = swap;
            }
        }
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}       
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("[ID %d][Iter %d] Comm Time: %f sec\n", id, num_iterations, time);

        num_iterations++;

        if( clock_gettime( CLOCK_REALTIME, &stop_iter) == -1 ) { perror( "clock gettime" );}       
        time_iter = (stop_iter.tv_sec - start_iter.tv_sec)+ (double)(stop_iter.tv_nsec - start_iter.tv_nsec)/1e9;
        printf("[ID %d][Iter %d] Iter Time: %f sec\n", id, num_iterations, time_iter);
    }
    // printf("I exit: %d %f %f\n", id, diff, diff_prev);

}
#endif

int main(){

#if  PREPROCESSING
    /* Preprocessing */
    string file_metadata = "/staging/vkp2/tye69227/256/graph.metadata";
    string file_matrix = "/staging/vkp2/tye69227/256/graph.txt";

    //vector<size_t> edgesDest;
    //vector<size_t> offsets;

    size_t *edgesDest, *offsets, *recip_offsets;

    size_t num_vertices;
    size_t num_edges;

    FILE *fp = fopen(file_metadata.c_str(), "r");
    fscanf(fp, "%lu %lu", &num_vertices, &num_edges);
    fclose(fp);

    edgesDest = (size_t *) calloc(num_edges , sizeof(size_t));
    offsets = (size_t *) calloc((num_vertices + 1) , sizeof(size_t));
    recip_offsets = (size_t *) calloc((num_vertices + 1) , sizeof(size_t));

//    edgesDest.resize(num_edges);
//    offsets.resize(num_vertices + 1);

    //offsets[num_vertices] = num_edges;

    read_matrix(file_matrix, edgesDest, offsets, recip_offsets, num_vertices + 1, num_edges);

/*
    cout << "offsets: \n";
    for (int k = 0; k <= num_vertices; k++)
        cout << offsets[k] << " ";
    cout << endl;

    cout << "recip_offsets: \n";
    for (int k = 0; k <= num_vertices; k++)
        cout << recip_offsets[k] << " ";
    cout << endl;

    cout << "edgesDest: \n";
    for (int k = 0; k < num_edges; k++)
        cout << edgesDest[k] << " ";
    cout << endl;
  */  


    read_offset_write(offsets,num_vertices);
    read_recipoffset_write(recip_offsets,num_vertices);

    read_partition_write(offsets, edgesDest, num_vertices);

    partition_meta(num_vertices, offsets);

    //read_partition(char filename[], size_t* edgesDest, size_t* offsets)

    // partition_meta(int edge_size, int offset_size)

//read_partition_write(size_t *offsets, size_t *edgesDest,int offset_size, int edge_size)

    //void read_partition_write_simple(size_t *edgesDest,int edge_size)
    // read_partition_write_simple(edgesDest,dump_count/NUM_WORKERS);

    /* End of Preprocessing*/

#else



    int tag = 50;
    MPI_Status status;  
    MPI_Init (NULL , NULL);
    int my_rank, proc_n;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);
    printf("%d\n",proc_n);

    //ring model
    ring_pagerank(my_rank,status,proc_n,tag);

    MPI_Finalize();


    // vector<size_t>::iterator ci;
    // for(int i=0; i<x; i++){
    //     for (ci = rows[i].begin(); ci != rows[i].end(); ci++)
    //         if(*ci != 0)
    //             printf("num: %f %d\n",*ci, i);
    // }

#endif

    return 0;
}
