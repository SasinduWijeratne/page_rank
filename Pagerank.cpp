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
#define PREPROCESSING 1
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

void read_partition_write(size_t *offsets, size_t *edgesDest,int offset_size, int edge_size){ // edge_size = one partition size
    FILE *f;
    char filename[20];
    int start_point = 0;
    int end_point = 0;
    int partition_number = 0;
    int global_count = 0;
    int internal_count = 0;

    while (end_point < offsets[offset_size] && partition_number < NUM_WORKERS){
        // printf("%d\n",partition_number);
        sprintf(filename, "partition%d.txt", partition_number);
        f = fopen(filename, "w");
        internal_count = 0;
        while((end_point - start_point) < edge_size){
            end_point = offsets[global_count];
            // printf("%d,%d,%d\n",start_point,end_point,edge_size); exit(0);
            if(end_point == offsets[offset_size]){
                // printf("Hit1\n");
                for(int kk= start_point; kk < end_point; kk++){
                    if(kk < edge_size*NUM_WORKERS)
                        fprintf(f, "%zu\n", edgesDest[kk]);
                }
                break;
            } 

            if((end_point - start_point) >= edge_size){
                if(internal_count == 0){
                    // printf("Hit2\n");
                    for(int kk= start_point; kk < end_point; kk++){
                        fprintf(f, "%zu\n", edgesDest[start_point + kk]);
                    }
                }
                else{
                    if((end_point - start_point) < 1.5*edge_size){
                        // printf("Hit3\n");
                        for(int kk= start_point; kk < end_point; kk++){
                            fprintf(f, "%zu\n", edgesDest[kk]);
                        }
                    }
                    else{
                        // printf("Hit4\n");
                        for(int kk= start_point; kk < offsets[global_count-1]; kk++){
                            fprintf(f, "%zu\n", edgesDest[kk]); 
                        }
                            global_count--;                           
                    }
                }
            }
                global_count++;
                internal_count++; 
        }
            
        partition_number++;
        start_point = end_point;
        fclose(f);
    }

}

void partition_meta(int edge_size, int offset_size, int per_vertex){
    FILE *f;
    char filename[20];
    sprintf(filename, "partition.metadata");
    f = fopen(filename, "w");
    fprintf(f, "%d %d %d\n", edge_size,offset_size,per_vertex);
    fclose(f);
}

void read_offset_write( size_t *offsets, int tot_size){ 
    FILE *f;
    char filename[20];
    sprintf(filename, "offset.txt");
    f = fopen(filename, "w");
    for(int i = 0; i < tot_size; i++ ){
        fprintf(f, "%zu\n", offsets[i]);
    }
    fclose(f);
}

void read_recipoffset_write( size_t *offsets, int tot_size){ 
    FILE *f;
    char filename[20];
    sprintf(filename, "recipoffset.txt");
    f = fopen(filename, "w");
    for(int i = 0; i < tot_size; i++ ){
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
        offsets[dest]++;
        if (prev != src) {
            prev++;
            for (; prev < src; prev++)
                recip_offsets[prev] = count;
            recip_offsets[src] = count;
            prev = src;
        }
        count++;
    }

    for(int kk=0; kk < num_vertices; kk++){
        if(kk > 0)
            offsets[kk] = offsets[kk-1] + offsets[kk];
            // printf("%d\n", offsets[kk]);
    }

    istream *infile2;
    infile2 = new ifstream(filename.c_str());
    string line2;

    while (getline(*infile2, line2)){
        stringstream lineStream(line2);
        size_t src, dest;
        lineStream >> src >> dest;
        if(dest > 0)
            edgesDest[offsets[dest-1] + temp_offsets[dest]] = src;
        else
            edgesDest[temp_offsets[dest]] = src;    
        temp_offsets[dest]++;
    }



    // for(int kk=count ; kk <num_edges; kk++ ){
    //    printf("%d\n",edgesDest[kk] );
    // }

    free(temp_offsets);

}
#if PREPROCESSING == 0
void ring_pagerank(int id, MPI_Status status, int proc_n, int tag,MPI_Request requests[]){

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


    string file_metadata = "partition.metadata";
    char filename[20];
    char filename1[20];
    char filename2[20];
    sprintf(filename, "partition%d.txt", id);
    sprintf(filename1, "offset.txt");
    sprintf(filename2, "recipoffset.txt");

    FILE *fp = fopen(file_metadata.c_str(), "r");
    fscanf(fp, "%lu %lu %lu", &num_edges,&tot_num_vertices, &num_vertices);
    fclose(fp);

    double *pr, *old_pr;
    size_t *edgesDest, *offsets, *recipoffsets;

    edgesDest = (size_t *) calloc(num_edges , sizeof(size_t));
    offsets = (size_t *) calloc((tot_num_vertices) , sizeof(size_t));
    recipoffsets = (size_t *) calloc((tot_num_vertices) , sizeof(size_t));

    pr = (double *) calloc(num_vertices , sizeof(double));
    old_pr = (double *) calloc(num_vertices , sizeof(double));

    read_partition_edges(filename,edgesDest);
    read_partition_offset(filename1,offsets);
    read_partition_offset(filename2,recipoffsets);

//    pr.resize(num_vertices);
//    old_pr.resize(num_vertices);

    pr[0] = 1;

    int num_iterations = 0;
    double diff = 1;
    double diff_prev = -1;
    int exit_flag = 0;
    // (fabs(diff - diff_prev) > convergence)

    // // iter_pagerank(vector< vector<size_t> > rows, vector<size_t> num_outgoing, vector<double> &pr, double &diff, int &num_iterations, int num_rows)
    while ((num_iterations < max_iterations) && (exit_flag == 0) ) {
        if(fabs(diff - diff_prev) <= convergence){
            // printf("Hit3 %d\n", id);
            tag = 51;
            exit_flag = 1;
        }
        diff_prev = diff;
        diff = 0;
        for(int sas_i = 0; sas_i < proc_n; sas_i++){
            // printf("tag: %d %d \n",id, tag);
            sum_pr = 0;
            dangling_pr = 0;

            for (size_t k = 0; k < num_vertices; k++) {
                double cpr = pr[k];
                sum_pr += cpr;
                if((num_vertices*((id + proc_n - sas_i)%proc_n)+k+1) < tot_num_vertices-1){
                    size_t temp_diff_offset = offsets[num_vertices*((id + proc_n - sas_i)%proc_n)+k+1] - offsets[num_vertices*((id + proc_n - sas_i)%proc_n)+k];
                    if (temp_diff_offset == 0) // Place 1
                        dangling_pr += cpr;
                }
           }


            diff += sum_pr;

            if (num_iterations == 0) {
                old_pr = pr;
            } else {
                for (size_t k = 0; k < num_vertices; k++) {
                    pr[k] /= sum_pr;
                }
                dangling_pr /= sum_pr;
                /* Normalize so that we start with sum equal to one */
                for (i = 0; i < num_vertices; i++) {
                    old_pr[i] = pr[i];
                }
            }
            /*
             * After normalisation the elements of the pagerank vector sum
             * to one
             */
            sum_pr = 1;
            
            /* An element of the A x I vector; all elements are identical */
            double one_Av = alpha * dangling_pr / (num_vertices*proc_n);

            /* An element of the 1 x I vector; all elements are identical */
            double one_Iv = (1 - alpha) * sum_pr / (num_vertices*proc_n);

            /* The difference to be checked for convergence */
            for (i = 0; i < num_vertices; i++) {
                /* The corresponding element of the H multiplication */
                double h = 0.0;
                int abs_i = num_vertices*((id + proc_n - sas_i)%proc_n) + i;
                if(abs_i < tot_num_vertices-1){
                    for (size_t ci = offsets[abs_i]; ci < offsets[abs_i+1]; ci++) {
                        /* The current element of the H vector */
                        if((edgesDest[ci] >= num_vertices*((id + proc_n - sas_i)%proc_n)) && (edgesDest[ci] < num_vertices*((id + proc_n - sas_i)%proc_n+1))){
                            double h_v = (recipoffsets[edgesDest[ci]+1]-recipoffsets[edgesDest[ci]])
                                        ? 1.0 / (recipoffsets[edgesDest[ci]+1]-recipoffsets[edgesDest[ci]])
                                        : 0.0;
                            h += h_v * old_pr[edgesDest[ci]%proc_n];                        
                        }
                    }

                }
                h *= alpha;
                pr[i] = h + one_Av + one_Iv;
            }

            //** if sas_i != proc_n 
            if (id%2 == 0) {
            
                // t1 = MPI_Wtime();  
                // printf("id: %d\n",id);
                // MPI_Send (&pr[0], num_vertices+1, MPI_INT, (id+1)%proc_n, tag, MPI_COMM_WORLD);
                // MPI_Recv (&pr[0], num_vertices+1, MPI_INT,(id == 0)? (proc_n-1): (id-1)%proc_n, tag, MPI_COMM_WORLD, &status);
                MPI_Isend(&pr[0], num_vertices, MPI_INT, (id+1)%proc_n, tag, MPI_COMM_WORLD, &requests[(id+1)%proc_n]);
                MPI_Wait(&requests[(id+1)%proc_n], &status);
                int use_id = (id == 0)? (proc_n-1): (id-1)%proc_n;
                MPI_Irecv(&pr[0], num_vertices, MPI_INT, use_id, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[use_id]);
                MPI_Wait(&requests[use_id], &status);
                if (status.MPI_TAG == 51){
                    // printf("Hit2 %d\n", id);
                    tag = 51;
                    exit_flag = 1;
                }
                
                // t2 = MPI_Wtime();
                // printf("\nRound trip(s): %f\n\n", t2-t1);    
            }
            else {
                // printf("id: %d\n",id);
                MPI_Irecv(&pr[0], num_vertices, MPI_INT,  (id-1)%proc_n, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[(id-1)%proc_n]);
                MPI_Wait(&requests[(id-1)%proc_n], &status);
                if (status.MPI_TAG == 51){
                    // printf("Hit2 %d\n", id);
                    tag = 51;
                    exit_flag = 1;
                }
                MPI_Isend(&pr[0], num_vertices, MPI_INT, (id+1)%proc_n, tag, MPI_COMM_WORLD, &requests[(id+1)%proc_n]);
                MPI_Wait(&requests[(id+1)%proc_n], &status);
            }

             MPI_Barrier(MPI_COMM_WORLD);

        }
            // printf("diff: %f\n",diff);
        // printf("%d %f\n",num_iterations,fabs(diff - diff_prev)/diff);
        num_iterations++;
        MPI_Barrier(MPI_COMM_WORLD);
        // if (trace) {
            // cout << num_iterations << ": ";
        // }
    }
    // printf("I exit: %d %f %f\n", id, diff, diff_prev);


}
#endif

int main(){

#if  PREPROCESSING
    /* Preprocessing */
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
    recip_offsets = (size_t *) calloc((num_vertices + 1) , sizeof(size_t));

//    edgesDest.resize(num_edges);
//    offsets.resize(num_vertices + 1);

    offsets[num_vertices] = num_edges;

    read_matrix(file_matrix, edgesDest, offsets, recip_offsets, num_vertices + 1, num_edges);

    read_offset_write(offsets,num_vertices);
    read_recipoffset_write(recip_offsets,num_vertices);

    read_partition_write(offsets, edgesDest,num_vertices, num_edges/NUM_WORKERS);

    partition_meta(num_edges, num_vertices+1,num_edges/NUM_WORKERS);

    //read_partition(char filename[], size_t* edgesDest, size_t* offsets)

    // partition_meta(int edge_size, int offset_size)

//read_partition_write(size_t *offsets, size_t *edgesDest,int offset_size, int edge_size)

    //void read_partition_write_simple(size_t *edgesDest,int edge_size)
    // read_partition_write_simple(edgesDest,dump_count/NUM_WORKERS);

    /* End of Preprocessing*/

#else



    int tag = 50;
    MPI_Status status;  
    MPI_Request request[NUM_WORKERS];
    MPI_Init (NULL , NULL);
    int my_rank, proc_n;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);
    // printf("%d\n",proc_n);

    //ring model
    ring_pagerank(my_rank,status,proc_n,tag, request);

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