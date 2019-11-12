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

// void iter_pagerank(vector< vector<size_t> > rows, vector<size_t> num_outgoing, vector<double> &pr, double &diff, int &num_iterations, int num_rows){

//     double sum_pr;
//     double dangling_pr;
//     double alpha = DEFAULT_ALPHA;
//     vector<double> old_pr;
//     vector<size_t>::iterator ci;
//     int i;

//     sum_pr = 0;
//     dangling_pr = 0;
    
//     for (size_t k = 0; k < pr.size(); k++) {
//         double cpr = pr[k];
//         sum_pr += cpr;
//         if (num_outgoing[k] == 0) {
//             dangling_pr += cpr;
//         }
//     }

//     if (num_iterations == 0) {
//         old_pr = pr;
//     } else {
//         /* Normalize so that we start with sum equal to one */
//         for (i = 0; i < pr.size(); i++) {
//             old_pr[i] = pr[i] / sum_pr;
//         }
//     }

//     /*
//      * After normalisation the elements of the pagerank vector sum
//      * to one
//      */
//     sum_pr = 1;
    
//     /* An element of the A x I vector; all elements are identical */
//     double one_Av = alpha * dangling_pr / num_rows;

//     /* An element of the 1 x I vector; all elements are identical */
//     double one_Iv = (1 - alpha) * sum_pr / num_rows;

//     /* The difference to be checked for convergence */
//     diff = 0;
//     for (i = 0; i < num_rows; i++) {
//         /* The corresponding element of the H multiplication */
//         double h = 0.0;
//         for (ci = rows[i].begin(); ci != rows[i].end(); ci++) {
//             /* The current element of the H vector */
//             double h_v = (num_outgoing[*ci])
//                 ? 1.0 / num_outgoing[*ci]
//                 : 0.0;
//             if (num_iterations == 0 && trace) {
//                 cout << "h[" << i << "," << *ci << "]=" << h_v << endl;
//             }
//             h += h_v * old_pr[*ci];
//         }
//         h *= alpha;
//         pr[i] = h + one_Av + one_Iv;
//         diff += fabs(pr[i] - old_pr[i]);
//     }
//     num_iterations++;
//     if (trace) {
//         cout << num_iterations << ": ";
//     }

// }

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

void row_pagerank(int id, MPI_Status status, int proc_n, int tag){
    vector<double> pr;

    double convergence = DEFAULT_CONVERGENCE;
    unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

    double sum_pr;
    double dangling_pr;
    double alpha = DEFAULT_ALPHA;
    vector<double> old_pr;
    vector<size_t>::iterator ci;
    int i;

    sum_pr = 0;
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
    old_pr.resize(len);
    pr[0] = 1;

    int num_iterations = 0;
    double diff = 1;

    // iter_pagerank(vector< vector<size_t> > rows, vector<size_t> num_outgoing, vector<double> &pr, double &diff, int &num_iterations, int num_rows)
    while (diff > convergence && num_iterations < max_iterations){
        // printf("id: %d\n",id);
        if(num_iterations > 0){

            if (id%2 == 0) {
            
                // t1 = MPI_Wtime();  
                // printf("id: %d\n",id);
                MPI_Send (&pr[0], len, MPI_INT, (id+1)%proc_n, tag, MPI_COMM_WORLD);
                MPI_Recv (&pr[0], len, MPI_INT,(id == 0)? (proc_n-1): (id-1)%proc_n, tag, MPI_COMM_WORLD, &status);
                
                // t2 = MPI_Wtime();
                // printf("\nRound trip(s): %f\n\n", t2-t1);    
            }
            else {
                // printf("id: %d\n",id);
                MPI_Recv (&pr[0], len, MPI_INT, (id-1)%proc_n, tag, MPI_COMM_WORLD, &status);
                MPI_Send (&pr[0], len, MPI_INT, (id+1)%proc_n, tag, MPI_COMM_WORLD); 
            }


        }
        if(pr[0] == -1) { 
            MPI_Send (&pr[0], len, MPI_INT, (id+1)%proc_n, tag, MPI_COMM_WORLD);

        break; }

        sum_pr = 0;
        dangling_pr = 0;

        for (size_t k = 0; k < pr.size(); k++) {
            double cpr = pr[k];
            sum_pr += cpr;
            if (num_outgoing[k] == 0) {
                dangling_pr += cpr;
            }
        }

        if (num_iterations == 0) {
            old_pr = pr;
        } else {
            /* Normalize so that we start with sum equal to one */
            for (i = 0; i < pr.size(); i++) {
                old_pr[i] = pr[i] / sum_pr;
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

        /* The difference to be checked for convergence */
        diff = 0;
        for (i = 0; i < num_rows; i++) {
            /* The corresponding element of the H multiplication */
            double h = 0.0;
            for (ci = rows[i].begin(); ci != rows[i].end(); ci++) {
                /* The current element of the H vector */
                double h_v = (num_outgoing[*ci])
                    ? 1.0 / num_outgoing[*ci]
                    : 0.0;
                // if (num_iterations == 0 && trace) {
                //     cout << "h[" << i << "," << *ci << "]=" << h_v << endl;
                // }
                h += h_v * old_pr[*ci];
            }
            h *= alpha;
            pr[i] = h + one_Av + one_Iv;
            diff += fabs(pr[i] - old_pr[i]);
            // printf("%f\n",dangling_pr);
            // if(num_iterations == 12){
            //     printf("h: %f one_Av: %f one_Iv: %f\n",h, one_Av, one_Iv);        
            // }
        }
        // printf("diff: %f\n",diff);
        num_iterations++;
        // if (trace) {
            // cout << num_iterations << ": ";
        // }
    }

    if (trace){
        double sum = 0;
        vector<double>::iterator cr;
        // cout << "(" << pr.size() << ") " << "[ ";
        for (cr = pr.begin(); cr != pr.end(); cr++) {
            // cout << *cr << " ";
            sum += *cr;
            // cout << "s = " << sum << " ";
        }
        cout << "] "<< sum << endl;
    }

    pr[0] = -1;
    if (id%2 == 0) {
    
        // t1 = MPI_Wtime();  

        MPI_Send (&pr[0], len, MPI_INT, (id+1)%proc_n, tag, MPI_COMM_WORLD);
        MPI_Recv (&pr[0], len, MPI_INT,(id == 0)? (proc_n-1): (id-1)%proc_n, tag, MPI_COMM_WORLD, &status);
        
        // t2 = MPI_Wtime();
        // printf("\nRound trip(s): %f\n\n", t2-t1);    
    }
    else {
    
        MPI_Recv (&pr[0], len, MPI_INT, (id-1)%proc_n, tag, MPI_COMM_WORLD, &status);
        MPI_Send (&pr[0], len, MPI_INT, (id+1)%proc_n, tag, MPI_COMM_WORLD);
    }


}


int main(){

    int tag = 50; /* Tag para as mensagens */
    MPI_Status status;  /* Status de retorno */

    MPI_Init (NULL , NULL);
    int my_rank, proc_n;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);
    // printf("%d\n",proc_n);

    row_pagerank(my_rank,status,proc_n,tag);

    MPI_Finalize();


    // vector<size_t>::iterator ci;
    // for(int i=0; i<x; i++){
    //     for (ci = rows[i].begin(); ci != rows[i].end(); ci++)
    //         if(*ci != 0)
    //             printf("num: %f %d\n",*ci, i);
    // }
    return 0;
}


