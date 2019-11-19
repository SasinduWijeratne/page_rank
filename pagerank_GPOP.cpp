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
/*mpic++ your_code_file.c
Execution

mpirun -np <no. of Processors> ./a.out
*/
#define PREPROCESSING 0

#define NUM_WORKERS 4

const double DEFAULT_ALPHA = 0.85;
const double DEFAULT_CONVERGENCE = 0.00001;
const unsigned long DEFAULT_MAX_ITERATIONS = 100;
int trace = 1;

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

void partition_meta(int edge_size, int offset_size, int per_vertex){
    FILE *f;
    char filename[25];
    sprintf(filename, "gpoppartition.metadata");
    f = fopen(filename, "w");
    fprintf(f, "%d %d %d\n", edge_size,offset_size,per_vertex);
    fclose(f);
}

// read_partition_write(offsets, edgesDest,num_edges, num_edges/NUM_WORKERS);
void read_partition_write(size_t *offsets, size_t *edgesDest,int offset_size, int edge_size, int num_vertices){ // edge_size = one partition size
    FILE *f;
    char filename[20];
    int start_point = 0;
    int end_point = 0;
    int partition_number = 0;
    int global_count = 0;
    int internal_count = 0;

    int per_worker = num_vertices/NUM_WORKERS;
    for(int i=0; i < NUM_WORKERS; i++){
        sprintf(filename, "gpoppartition%d.txt", i);
        f = fopen(filename, "w");
        int start = offsets[i*per_worker];
        int end =  offsets[(i+1)*per_worker];
        for(int j=start; j < end; j++){
            fprintf(f, "%zu\n", edgesDest[j]);
        }
        fclose(f);
    }

}

void read_offset_write( size_t *offsets, int tot_size){ 
    FILE *f;
    char filename[20];
    sprintf(filename, "gpopoffset.txt");
    f = fopen(filename, "w");
    for(int i = 0; i < tot_size; i++ ){
        fprintf(f, "%zu\n", offsets[i]);
    }
    fclose(f);
}

void read_vector(const string &filename, int &len, vector<size_t> &num_outgoing){

    istream *infile;

    infile = new ifstream(filename.c_str());

    string line;

    getline(*infile, line);
    len = strtol(line.c_str(), NULL, 10);
    num_outgoing.resize(len);
    int count = 0;
    while (getline(*infile, line)){
        int num = strtol(line.c_str(), NULL, 10);
        num_outgoing[count] = num;
        count++;
    }
}

void read_matrix(const string &filename, size_t *edgesDest, size_t *offsets){
    istream *infile;
    infile = new ifstream(filename.c_str());
    string line;

    size_t count = 0;
    size_t prev = -1;
    while (getline(*infile, line)){
        stringstream lineStream(line);
        size_t src, dest;
        lineStream >> src >> dest;
        if (prev != src) {
            prev++;
            for (; prev < src; prev++)
                offsets[prev] = count;
            offsets[src] = count;
            prev = src;
        }
        edgesDest[count] = dest;
        count++;
    }
}
#if PREPROCESSING == 0
void gpop_pagerank(int id, MPI_Status status, int proc_n, int tag){

    double convergence = DEFAULT_CONVERGENCE;
    unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

    double sum_pr;
    double dangling_pr;
    double alpha = DEFAULT_ALPHA;

    int i,j;

    sum_pr = 0;
    dangling_pr = 0;

    size_t num_vertices;
    size_t tot_num_vertices;
    size_t num_edges;

    // string file_vector = "vector" + to_string(id) +".txt";
    // string file_matrix = "mat" + to_string(id) +".txt";
    string file_metadata = "gpoppartition.metadata";
    char filename[20];
    char filename1[20];
    sprintf(filename, "gpoppartition%d.txt", id);
    sprintf(filename1,"gpopoffset.txt");

    FILE *fp = fopen(file_metadata.c_str(), "r");
    fscanf(fp, "%lu %lu %lu", &num_edges,&tot_num_vertices, &num_vertices);
    fclose(fp);

    double  *pr, *old_pr;
    size_t *edgesDest, *offsets;

    edgesDest = (size_t *) calloc(num_edges , sizeof(size_t));
    offsets = (size_t *) calloc((tot_num_vertices) , sizeof(size_t));

    pr = (double *) calloc(num_vertices , sizeof(double));
    old_pr = (double *) calloc(num_vertices , sizeof(double));

    read_partition_edges(filename,edgesDest);
    read_partition_offset(filename1,offsets);

    int num_iterations = 0;
    double diff = 99;
    double diff_i = 1;
    double diff_i_old = -1;
    old_pr[0] = 1;
    pr[0] = 1;

    double *sending_arr = (double *)malloc((tot_num_vertices + NUM_WORKERS) * sizeof(double)); 
    double *receiving_arr = (double *)malloc((tot_num_vertices + NUM_WORKERS) * sizeof(double)); 


    while (fabs(diff_i) > convergence && num_iterations < max_iterations && fabs(diff_i_old - diff_i) > 0.0000001){
        if(num_iterations == 0){
            double temp_init =   1.0/tot_num_vertices;
            if(temp_init > 0.00000000000000000001){
                for (i = 0; i < num_vertices; ++i)
                    old_pr[i] = temp_init;  
            }
            else{
                old_pr[0] = 1;
            }

        }
    
    // Scatter Phase

        for(j=0; j<tot_num_vertices; j++){
            sending_arr[j] = 0;
        }
        for (i = 0; i < num_vertices; i++) {
            double temp_pr = (offsets[i+1]-offsets[i]) ? 1.0 / (offsets[i+1]-offsets[i]) : 0.0;

            for (size_t ci = offsets[i]; ci < offsets[i+1]; ci++) {
                sending_arr[edgesDest[ci]] += (temp_pr*old_pr[i]);
            }
        }

    for(i = 0; i < NUM_WORKERS; i++)
        sending_arr[num_vertices*(i+1)+i] = diff;

    int status_all2all = MPI_Alltoall(&sending_arr[0], num_vertices+1, MPI_DOUBLE, &receiving_arr[0], num_vertices+1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    diff_i_old = diff_i;
    diff_i = 0;
    // printf("id: %d diff %f\n",id, diff );
    for(i = 0; i < NUM_WORKERS; i++){
        diff_i += receiving_arr[(i+1)*num_vertices+i];
    }
    // printf("id: %d i_diff: %f\n",id, diff_i);

    //Gather Phase
    sum_pr = 0;
    for (size_t k = 0; k < num_vertices; k++) {
        double cpr = pr[k];
        sum_pr += cpr;
    }
    if(sum_pr == 0) sum_pr = 1; // Hack for test data
    diff = 0;
    for(i=0; i < num_vertices; i++){
        pr[i] = 0;
        for(j=0; j < NUM_WORKERS; j++){
            pr[i] += receiving_arr[j *(num_vertices+1) +i];
        }
            pr[i] = ((1-alpha)/tot_num_vertices) + alpha*pr[i];

        diff += fabs(pr[i] - old_pr[i]/sum_pr);
    }

    for(i = 0; i < num_vertices; i++)
        old_pr[i] = pr[i];

    num_iterations++;

    }
    // for ( i = 0; i < pr.size(); ++i)
    // {
    //     printf("%f\n",pr[i]);
    // }

}
#endif

int main(){
#if PREPROCESSING

    string file_metadata = "graph.metadata";
    string file_matrix = "graph.txt";

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

    read_matrix(file_matrix,edgesDest,offsets);
    read_offset_write(offsets,num_vertices);
    read_partition_write(offsets, edgesDest,num_edges, num_edges/NUM_WORKERS, num_vertices);
    partition_meta(num_edges, num_vertices,num_vertices/NUM_WORKERS);

#else
    int tag = 50;
    MPI_Status status;  

    MPI_Init (NULL , NULL);
    int my_rank, proc_n;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);
    // printf("%d\n",proc_n);

    //ring model
    gpop_pagerank(my_rank,status,proc_n,tag);

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



