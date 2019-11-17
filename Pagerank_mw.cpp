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

const double DEFAULT_ALPHA = 0.85;
const double DEFAULT_CONVERGENCE = 0.00001;
const unsigned long DEFAULT_MAX_ITERATIONS = 10000;
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

void read_matrix(const string &filename, size_t *edgesDest, size_t *offsets, int num_vertices, int num_edges){
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
    MPI_Request requests[16];
    int flag;

    MPI_Init (&argc , &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);
    p = proc_n - 1;


    string file_metadata = "partition.metadata";

    size_t num_vertices;
    size_t tot_num_vertices;
    size_t num_edges;

    FILE *fp = fopen(file_metadata.c_str(), "r");
    fscanf(fp, "%lu %lu %lu", &num_edges,&tot_num_vertices, &num_vertices);
    fclose(fp);

    double *pr, *old_pr;
    size_t *edgesDest, *offsets;

    pr = (double *) calloc(tot_num_vertices , sizeof(double));
    old_pr = (double *) calloc(tot_num_vertices , sizeof(double));

    if(my_rank > 0){
        char filename1[20];
        sprintf(filename1, "offset.txt");
        edgesDest = (size_t *) malloc(num_edges * sizeof(size_t));
        offsets = (size_t *) malloc((tot_num_vertices) * sizeof(size_t));
        char filename[20];
        sprintf(filename, "partition%d.txt", my_rank-1);
        read_partition_edges(filename,edgesDest);
        read_partition_offset(filename1,offsets);
        

    }

    double convergence = DEFAULT_CONVERGENCE;
    unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

    double sum_pr;
    double dangling_pr;
    double alpha = DEFAULT_ALPHA;

    sum_pr = 0;
    dangling_pr = 0;
    pr[0] = 1;

    int num_iterations = 0;
    double diff = 1;
    double diff_prev = -1;
    int len = num_vertices;
    int num_rows = tot_num_vertices;

    /********************************/

    MPI_Barrier(MPI_COMM_WORLD);

    while ((fabs(diff - diff_prev) > convergence) && num_iterations < max_iterations){
        if (my_rank == 0) {
            double t1, t2;
            t1 = MPI_Wtime();
            for (i = 1; i < proc_n; i++)
                MPI_Isend(&pr[0], tot_num_vertices, MPI_INT, i, tag, MPI_COMM_WORLD, &requests[i]);
            for (i = 1; i < proc_n; i++)
                MPI_Wait(&requests[i], &status);
            t2 = MPI_Wtime();
            printf("[Master] Time for MPI recv: %f\n", t2-t1);
        } else {
            double t1, t2;
            t1 = MPI_Wtime();
            MPI_Irecv(&pr[0], tot_num_vertices, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            t2 = MPI_Wtime();
            printf("[Node %d] Time for MPI recv: %f\n", my_rank, t2-t1);
        }

        if(my_rank > 0) {
            sum_pr = 0;
            dangling_pr = 0;

            for (size_t k = 0; k < num_vertices; k++) {
                double cpr = pr[k];
                sum_pr += cpr;
                if (offsets[k+1] - offsets[k] == 0) {
                    dangling_pr += cpr;
                }
            }
            diff_prev = diff;
            if (num_iterations == 0) {
                old_pr[0] = 1;
                diff = 999;
            } else {
                /* Normalize so that we start with sum equal to one */
                diff = 0;
                
                double tmp;
                for (i = 0; i < num_vertices; i++) {
                    tmp = pr[i] / sum_pr;
                    diff += fabs(tmp - old_pr[i]);
                    old_pr[i] = tmp;
                }
            }
            /*
             * After normalisation the elements of the pagerank vector sum
             * to one
             */
            sum_pr = 1;
            
             /* An element of the A x I vector; all elements are identical */
            double one_Av = alpha * dangling_pr / num_rows;

            /* An element of the 1 x I vector; all elements are identical */
            double one_Iv = (1 - alpha) * sum_pr / num_rows;

            for (i = 0; i < tot_num_vertices; i++) {
                /* The corresponding element of the H multiplication */
                double h = 0.0;
                for (size_t ci = offsets[i]; ci < offsets[i+1]; ci++) {
                    /* The current element of the H vector */
                    double h_v = (offsets[edgesDest[ci]+1]-offsets[edgesDest[ci]])
                        ? 1.0 / (offsets[edgesDest[ci]+1]-offsets[edgesDest[ci]])
                        : 0.0;
                    h += h_v * old_pr[edgesDest[ci]];
                }
                h *= alpha;
                pr[i] = h + one_Av + one_Iv;
            }

            
        }

        if (my_rank == 0) {
            double t1, t2;
            t1 = MPI_Wtime();
            for (i = 0; i < p; i++)
                MPI_Irecv(&pr[i*len], len, MPI_INT, i+1, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[i]);
            for (i = 0; i < p; i++)
                MPI_Wait(&requests[i], &status);
            t2 = MPI_Wtime();
            if (status.MPI_TAG == 51) {
                diff = 0;
                diff_prev = 0;
            }
            printf("[Master] Time for MPI recv: %f\n", t2-t1);
        } else {
            double t1, t2;
            if (fabs(diff - diff_prev) <= convergence)
                tag = 51;
            t1 = MPI_Wtime();
            MPI_Isend(&pr[(my_rank-1)*len], len, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            t2 = MPI_Wtime();
            printf("[Node %d] Time for MPI send: %f\n", my_rank, t2-t1);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        num_iterations++;
    }


    MPI_Finalize();
    return 0;
}

