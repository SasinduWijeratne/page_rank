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

#define NUM_WORKERS 16

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
        sprintf(filename, "gpoppartition%d.txt", pid);
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
    char filename[30];

    for (int pid = 0; pid < NUM_WORKERS; pid++) {
        sprintf(filename, "gpoppartition%d.metadata", pid);
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
        sprintf(filename, "gpopoffset%d.txt", pid);
        f = fopen(filename, "w");
        size_t left, right;
        left = num_vertices_per_worker*pid;
        right = (pid == NUM_WORKERS - 1) ? tot_size : num_vertices_per_worker*(pid+1);
        for (int i = left; i <= right; i++)
            fprintf(f, "%zu\n", offsets[i] - offsets[left]);  
        fclose(f);
    }
}

void read_matrix(const string &filename, size_t *edgesDest, size_t *offsets, size_t size){
    istream *infile;
    infile = new ifstream(filename.c_str());
    string line;

    size_t count = 0;
    size_t prev = -1;
    size_t src, dest;
    while (getline(*infile, line)){
        stringstream lineStream(line);
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
/*
    for(size_t i = 0; i < size; i++)
        if(offsets[i] == 0)
            offsets[i] = offsets[i-1];
*/
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
    char file_metadata[30];
    sprintf(file_metadata, "gpoppartition%d.metadata", id);
    char filename[30];
    char filename1[30];
    sprintf(filename, "gpoppartition%d.txt", id);
    sprintf(filename1,"gpopoffset%d.txt",id);

    cout << file_metadata << " " << filename << " " << filename1 << endl;

    FILE *fp = fopen(file_metadata, "r");
    fscanf(fp, "%lu %lu %lu", &tot_num_vertices, &num_edges, &num_vertices);
    fclose(fp);

    cout << tot_num_vertices << " " << num_edges << " " << num_vertices << " " << id << endl;

    double  *pr, *old_pr;
    size_t *edgesDest, *offsets;

    edgesDest = (size_t *) calloc(num_edges , sizeof(size_t));
    offsets = (size_t *) calloc(num_vertices , sizeof(size_t));

    pr = (double *) calloc(num_vertices , sizeof(double));
    old_pr = (double *) calloc(num_vertices , sizeof(double));

    cout << "Alloc Success " << id << "\n";

    read_partition_edges(filename,edgesDest);
    read_partition_offset(filename1,offsets);

    cout << "Reading Success " << id << "\n";

    int num_iterations = 0;
    double diff = 99;
    double diff_i = 1;
    double diff_i_old = -1;
    old_pr[0] = 1;
    pr[0] = 1;

    size_t denp = (tot_num_vertices/NUM_WORKERS);
    size_t per_buffer_size = tot_num_vertices - (NUM_WORKERS-1)*denp;
    size_t max_buffer_size = per_buffer_size*NUM_WORKERS + NUM_WORKERS;

    // printf("per_buffer_size = %zu %zu\n",per_buffer_size,max_buffer_size);

    cout << denp << " " << per_buffer_size << " " << max_buffer_size << " " << id << endl;

    double *sending_arr, *receiving_arr; 

    cout << "malloc size " << max_buffer_size * sizeof(double) << " " << id << endl;
    sending_arr = (double *)malloc(max_buffer_size * sizeof(double)); 
    cout << "hi " << id << "\n";
    if (sending_arr == NULL) cout << "sending_arr Error " << id << "\n";

    receiving_arr = (double *)malloc(max_buffer_size * sizeof(double)); 
    cout << "hi2 " << id << "\n";
    if (receiving_arr == NULL) cout << "recving_arr Error " << id << "\n";

    cout << "Alloc2 Success " << id << "\n";

    struct timespec start_iter, stop_iter, start, stop;
    double time_iter, time;

    while (num_iterations < max_iterations && fabs(diff_i_old - diff_i) > 0.00001){

        if(num_iterations == 0){
            double temp_init = 1.0/tot_num_vertices;
            if(temp_init > 0.00000000000000000001){
                for (i = 0; i < num_vertices; ++i)
                    old_pr[i] = temp_init;  
            }
            else{
                old_pr[0] = 1;
            }

        }

    // Scatter Phase

        clock_gettime(CLOCK_REALTIME, &start_iter);

        clock_gettime(CLOCK_REALTIME, &start);

        for(j=0; j<max_buffer_size; j++){
            sending_arr[j] = 0;
        }
        
        for (i = 0; i < num_vertices; i++) {
            double temp_pr = (offsets[i+1]-offsets[i]) ? 1.0 / (offsets[i+1]-offsets[i]) : 0.0;

            for (size_t ci = offsets[i]; ci < offsets[i+1]; ci++) {
                size_t desination = edgesDest[ci];
                size_t partition = edgesDest[ci]/denp;
                size_t relative = (partition == NUM_WORKERS-1) ? (partition*per_buffer_size + (desination - (partition)*denp)): (partition*per_buffer_size +desination%denp);
                sending_arr[relative] += (temp_pr*old_pr[i]);
            }
        }

        for(size_t send_id = 0; send_id < NUM_WORKERS; send_id++){
            sending_arr[(per_buffer_size)*(send_id+1) + send_id] = diff;
        }
        
        clock_gettime(CLOCK_REALTIME, &stop);
        time = (stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("[ID %d][Iter %d] Comp1 Time: %f sec\n", id, num_iterations, time);

        clock_gettime(CLOCK_REALTIME, &start);
        
        int status_all2all = MPI_Alltoall(sending_arr, per_buffer_size+1, MPI_DOUBLE, receiving_arr, per_buffer_size+1, MPI_DOUBLE, MPI_COMM_WORLD);

        clock_gettime(CLOCK_REALTIME, &stop);
        time = (stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("[ID %d][Iter %d] Comm Time: %f sec\n", id, num_iterations, time);
        
        clock_gettime(CLOCK_REALTIME, &start);

        diff_i_old = diff_i;
        diff_i = 0;
   
        for(size_t send_id = 0; send_id < NUM_WORKERS; send_id++){
            diff_i += receiving_arr[(per_buffer_size)*(send_id+1) + send_id];
        }
           
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
                pr[i] += receiving_arr[j *(per_buffer_size+1) +i];
            }
            pr[i] = ((1-alpha)/tot_num_vertices) + alpha*pr[i];

            diff += fabs(pr[i] - old_pr[i]);
        }

        for(i = 0; i < num_vertices; i++)
            old_pr[i] = pr[i];
        
        clock_gettime(CLOCK_REALTIME, &stop);
        time = (stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("[ID %d][Iter %d] Comp2 Time: %f sec\n", id, num_iterations, time);
        
        clock_gettime(CLOCK_REALTIME, &stop_iter);
        time_iter = (stop_iter.tv_sec - start_iter.tv_sec) + (double)(stop_iter.tv_nsec - start_iter.tv_nsec)/1e9;
        printf("[ID %d][Iter %d] Total Time: %f sec\n", id, num_iterations, time_iter);
        printf("[ID %d][Iter %d] Diff: %f \n", id, num_iterations, diff);
        
        num_iterations++;
    }

}
#endif

int main(int argc, char** argv){
#if PREPROCESSING

    string file_metadata = "/staging/vkp2/tye69227/256/graph.metadata";
    string file_matrix = "/staging/vkp2/tye69227/256/graph.txt";

    size_t *edgesDest, *offsets, *recip_offsets;

    size_t num_vertices;
    size_t num_edges;

    FILE *fp = fopen(file_metadata.c_str(), "r");
    fscanf(fp, "%lu %lu", &num_vertices, &num_edges);
    fclose(fp);

    cout << num_vertices << " " << num_edges << endl;

    edgesDest = (size_t *) calloc(num_edges , sizeof(size_t));
    offsets = (size_t *) calloc((num_vertices + 1) , sizeof(size_t));

    read_matrix(file_matrix,edgesDest,offsets, num_vertices + 1);
    offsets[num_vertices] = num_edges;

/*    
    cout << "offsets:\n";
    for (size_t k = 0; k <= num_vertices; k++)
        cout << offsets[k] << " ";
    cout << endl;
    
    cout << "edgesDest:\n";
    for (size_t k = 0; k < num_edges; k++)
        cout << edgesDest[k] << " ";
    cout << endl;
*/

    read_offset_write(offsets,num_vertices);
    read_partition_write(offsets, edgesDest, num_vertices);
    partition_meta(num_vertices,offsets);

#else
    int tag = 50;
    MPI_Status status;  
    
    MPI_Init (&argc , &argv);
    int my_rank, proc_n;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);

    gpop_pagerank(my_rank,status,proc_n,tag);

    MPI_Finalize();

#endif

    return 0;
}



