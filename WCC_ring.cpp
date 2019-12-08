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
    size_t num_edges, num_edges2;

    MPI_Request request;

    char file_metadata[30];
    char file_metadata2[30];
    char filename[30];
    char filename1[30];
    char filename2[30];
    char filename3[30];

    sprintf(file_metadata, "partition%d.metadata", id);
    sprintf(file_metadata2, "gpoppartition%d.metadata", id);
    sprintf(filename, "partition%d.txt", id);
    sprintf(filename1, "offset%d.txt", id);
    sprintf(filename2, "gpoppartition%d.txt", id);
    sprintf(filename3, "gpopoffset%d.txt", id);

    FILE *fp = fopen(file_metadata, "r");
    fscanf(fp, "%lu %lu %lu", &tot_num_vertices, &num_edges2, &num_vertices);
    fclose(fp);

    fp = fopen(file_metadata2, "r");
    fscanf(fp, "%lu %lu %lu", &tot_num_vertices, &num_edges, &num_vertices);
    fclose(fp);

    size_t *pr, *old_pr, *rank;
    size_t *edgesDest, *offsets, *edgesDest2, *offsets2;
    size_t left, right;

    edgesDest = (size_t *) calloc(num_edges , sizeof(size_t));
    offsets = (size_t *) calloc((tot_num_vertices) , sizeof(size_t));
    edgesDest2 = (size_t *) calloc(num_edges2 , sizeof(size_t));
    offsets2 = (size_t *) calloc((tot_num_vertices) , sizeof(size_t));

    pr = (size_t *) calloc(tot_num_vertices/NUM_WORKERS + NUM_WORKERS, sizeof(size_t));
    old_pr = (size_t *) calloc(tot_num_vertices, sizeof(size_t));
    rank = (size_t *) calloc(tot_num_vertices/NUM_WORKERS + NUM_WORKERS, sizeof(size_t));

    size_t *pr1;
    pr1 = (size_t *) malloc((tot_num_vertices/NUM_WORKERS + NUM_WORKERS) * sizeof(size_t));

    read_partition_edges(filename,edgesDest2);
    read_partition_offset(filename1,offsets2);
    read_partition_edges(filename2,edgesDest);
    read_partition_offset(filename3,offsets);

    left = (tot_num_vertices/NUM_WORKERS)*id;

    for (size_t v = 0; v < num_vertices; v++)
        pr[v] = tot_num_vertices;
    for (size_t v = 0; v < tot_num_vertices; v++)
        old_pr[v] = tot_num_vertices;
    for (size_t v = 0; v < num_vertices; v++)
        rank[v] = left+v;

    int num_iterations = 0;

    // (fabs(diff - diff_prev) > convergence)
    struct timespec start_iter, stop_iter, start, stop;
    double time_iter, time;
    int hasUpdate = 1;

    while ((num_iterations < max_iterations) && hasUpdate) {

        int count = 0;
        for (size_t v = 0; v < tot_num_vertices; v++)
            if (old_pr[v] != 0)
                count++;
        cout << "Non-zero: " << count << endl;


        if( clock_gettime( CLOCK_REALTIME, &start_iter) == -1 ) { perror( "clock gettime" );}        
        if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
        

        int update_count = 0;
        for (size_t v = 0; v < num_vertices; v++) {
            for (size_t ci = offsets2[v]; ci < offsets2[v+1]; ci++) {
                size_t u = edgesDest2[ci];
                if (old_pr[u] < rank[v]) 
                    rank[v] = old_pr[u];
            }
            if (pr[v] != rank[v]) {
                update_count++;
            }
            pr[v] = rank[v];
        }

        hasUpdate = 0;
        if (update_count > 0)
            hasUpdate = 1;

        pr[tot_num_vertices/NUM_WORKERS + NUM_WORKERS - 1] = hasUpdate;

        hasUpdate = 0;
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}       
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("[ID %d][Round 1][Iter %d] Comp Time: %f sec\n", id, num_iterations, time);
        
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
                hasUpdate |= pr[tot_num_vertices/NUM_WORKERS + NUM_WORKERS - 1];

                MPI_Isend(&pr[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_UNSIGNED_LONG, (id+1)%proc_n, tag, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                int use_id = (id == 0)? (proc_n-1): (id-1)%proc_n;

                MPI_Irecv(&pr[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_UNSIGNED_LONG, use_id, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
            }
            else {

                MPI_Irecv(&pr1[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_UNSIGNED_LONG,  (id-1)%proc_n, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);

                for (size_t k = 0; k < current_partition_size; k++)
                    old_pr[left + k] = pr[k];
                hasUpdate |= pr[tot_num_vertices/NUM_WORKERS + NUM_WORKERS - 1];

                MPI_Isend(&pr[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_UNSIGNED_LONG, (id+1)%proc_n, tag, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);

                size_t * swap = pr;
                pr = pr1;
                pr1 = swap;
            }
        }
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}       
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("[ID %d][Round 1][Iter %d] Comm Time: %f sec\n", id, num_iterations, time);

        if( clock_gettime( CLOCK_REALTIME, &stop_iter) == -1 ) { perror( "clock gettime" );}       
        time_iter = (stop_iter.tv_sec - start_iter.tv_sec)+ (double)(stop_iter.tv_nsec - start_iter.tv_nsec)/1e9;
        printf("[ID %d][Round 1][Iter %d] Iter Time: %f sec\n", id, num_iterations, time_iter);
        num_iterations++;
    }

    num_iterations = 0;
    hasUpdate = 1;
    while ((num_iterations < max_iterations) && (hasUpdate || num_iterations < 3)) {
        if( clock_gettime( CLOCK_REALTIME, &start_iter) == -1 ) { perror( "clock gettime" );}        
        if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
        
        hasUpdate = 0;
        for (size_t v = 0; v < num_vertices; v++) {
            for (size_t ci = offsets[v]; ci < offsets[v+1]; ci++) {
                size_t u = edgesDest[ci];
                if (old_pr[u] < rank[v]) 
                    rank[v] = old_pr[u];
            }
            if (pr[v] != rank[v]) 
                hasUpdate = 1;
            pr[v] = rank[v];
        }
        
        pr[tot_num_vertices/NUM_WORKERS + NUM_WORKERS - 1] = hasUpdate;

        hasUpdate = 0;
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}       
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("[ID %d][Round 2][Iter %d] Comp Time: %f sec\n", id, num_iterations, time);
        
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
                hasUpdate |= pr[tot_num_vertices/NUM_WORKERS + NUM_WORKERS - 1];

                MPI_Isend(&pr[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_UNSIGNED_LONG, (id+1)%proc_n, tag, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                int use_id = (id == 0)? (proc_n-1): (id-1)%proc_n;

                MPI_Irecv(&pr[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_UNSIGNED_LONG, use_id, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
            }
            else {

                MPI_Irecv(&pr1[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_UNSIGNED_LONG,  (id-1)%proc_n, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);

                for (size_t k = 0; k < current_partition_size; k++)
                    old_pr[left + k] = pr[k];
                hasUpdate |= pr[tot_num_vertices/NUM_WORKERS + NUM_WORKERS - 1];

                MPI_Isend(&pr[0], tot_num_vertices/NUM_WORKERS + NUM_WORKERS, MPI_UNSIGNED_LONG, (id+1)%proc_n, tag, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);

                size_t *swap = pr;
                pr = pr1;
                pr1 = swap;
            }
        }
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}       
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("[ID %d][Round 2][Iter %d] Comm Time: %f sec\n", id, num_iterations, time);

        if( clock_gettime( CLOCK_REALTIME, &stop_iter) == -1 ) { perror( "clock gettime" );}       
        time_iter = (stop_iter.tv_sec - start_iter.tv_sec)+ (double)(stop_iter.tv_nsec - start_iter.tv_nsec)/1e9;
        printf("[ID %d][Round 2][Iter %d] Iter Time: %f sec\n", id, num_iterations, time_iter);
        num_iterations++;
    }


}

int main(){
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

    return 0;
}
