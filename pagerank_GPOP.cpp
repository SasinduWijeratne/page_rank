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
#define NUM_WORKERS 4
#define _PER_VERTICES 7499
#define _PER_VERTICESX2 _PER_VERTICES*2
#define NUM_VERTICES _PER_VERTICES*NUM_WORKERS

const double DEFAULT_ALPHA = 0.85;
const double DEFAULT_CONVERGENCE = 0.00001;
const unsigned long DEFAULT_MAX_ITERATIONS = 10000;
int trace = 1;

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

void read_matrix(const string &filename, int len, vector< vector<size_t> > &rows){

    istream *infile;

    infile = new ifstream(filename.c_str());

    string line;

    rows.resize(len);
    int count = 0;
    while (getline(*infile, line)){
        stringstream  lineStream(line);
        size_t value;

        while(lineStream >> value)
        {
            rows[count].push_back(value);
            // if(value != 0)
            //     printf("%d\n",value);

        }
        count++;
    }
}

void gpop_pagerank(int id, MPI_Status status, int proc_n, int tag){
    vector<double> pr;
    vector<double> count_pr;

    double convergence = DEFAULT_CONVERGENCE;
    unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

    double sum_pr;
    double dangling_pr;
    double alpha = DEFAULT_ALPHA;
    vector<double> old_pr;
    vector<size_t>::iterator ci;
    int i,j;

    sum_pr = 1;
    dangling_pr = 0;

    // string file_vector = "vector" + to_string(id) +".txt";
    // string file_matrix = "mat" + to_string(id) +".txt";
    string file_vector = "vector0.txt";
    string file_matrix = "mat0.txt";

    vector< vector<size_t> > rows;
    vector<size_t> num_outgoing;
    int len;

    read_vector(file_vector,len, num_outgoing);
    read_matrix(file_matrix,len,rows);

    int num_rows = len;
    pr.resize(len);
    count_pr.resize(len);
    old_pr.resize(len);
    int num_iterations = 0;
    double diff = 1;
    double diff_i = 1;
    old_pr[0] = 1;
    pr[0] = 1;

    double *sending_arr = (double *)malloc((NUM_WORKERS * _PER_VERTICESX2 + NUM_WORKERS) * sizeof(double *)); 
    double *receiving_arr = (double *)malloc((NUM_WORKERS * _PER_VERTICESX2 + NUM_WORKERS) * sizeof(double *)); 


    while (diff_i > convergence && num_iterations < max_iterations){

        if(num_iterations == 0){
            for (i; i < old_pr.size(); ++i)
            {
                old_pr[i] = 1/num_outgoing[i];
            }           
        }
    
    // Scatter Phase
    for(i=0; i<NUM_WORKERS; i++){
        for(j=0; j<_PER_VERTICESX2; j++){
            sending_arr[i * _PER_VERTICESX2 + j] = 0;
        }
    }
    for (i = 0; i < _PER_VERTICES; i++) {
        /* The corresponding element of the H multiplication */
        double deg = num_outgoing[i];
        for (ci = rows[i].begin(); ci != rows[i].end(); ci++) {
            int pos_ = *ci;
            double temp_pr = (deg) ? 1.0 / deg : 0.0;
            // if (num_iterations == 0 && trace) {
            //     cout << "h[" << i << "," << *ci << "]=" << h_v << endl;
            // }
            int scatter_partition = pos_/_PER_VERTICES;
            int relative_pos = pos_%_PER_VERTICES;
            sending_arr[scatter_partition * _PER_VERTICESX2 + 2*relative_pos] += (temp_pr*(old_pr[i]/sum_pr));

            sending_arr[scatter_partition * _PER_VERTICESX2 + 2*relative_pos+1] += 1;
        }

    }
    sending_arr[_PER_VERTICESX2*NUM_WORKERS + 1] = diff;
    MPI_Barrier(MPI_COMM_WORLD);
    int status_all2all = MPI_Alltoall(sending_arr, _PER_VERTICESX2+1, MPI_DOUBLE, receiving_arr, _PER_VERTICESX2+1, MPI_DOUBLE, MPI_COMM_WORLD);
    diff_i = 0;
    for(i = 0; i < NUM_WORKERS; i++){
        diff_i += receiving_arr[i*_PER_VERTICESX2 + 1];
    }

    //Gather Phase
    old_pr = pr;
    sum_pr = 0;
    for (size_t k = 0; k < pr.size(); k++) {
        double cpr = pr[k];
        sum_pr += cpr;
    }
    if(sum_pr == 0) sum_pr = 1; 
    diff = 0;
    for(i=0; i < _PER_VERTICES; i++){
        pr[i] = 0;
        count_pr[i] = 0;
        for(j=0; j < NUM_WORKERS; j++){
            pr[i] += receiving_arr[j *_PER_VERTICESX2 +2*i];
            count_pr[i] += receiving_arr[j *_PER_VERTICESX2 +2*i+1];
        }
        if(count_pr[i] > 0)
            pr[i] = ((1-alpha)/count_pr[i]) + alpha*pr[i];
        else
            pr[i] = alpha*pr[i];
        diff += fabs(pr[i] - old_pr[i]/sum_pr);
    }
    num_iterations++;
    MPI_Barrier(MPI_COMM_WORLD);


    }
    // for ( i = 0; i < pr.size(); ++i)
    // {
    //     printf("%f\n",pr[i]);
    // }

}

int main(){

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
    return 0;
}



