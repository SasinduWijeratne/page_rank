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
    char filename[20];
    int *message;

    MPI_Status status;  /* Status de retorno */
    MPI_Request request;
    MPI_Request requests[16];
    int flag;

    MPI_Init (&argc , &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);
    p = proc_n - 1;


    /*local var */
    vector<double> pr;

    double convergence = DEFAULT_CONVERGENCE;
    unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

    double sum_pr;
    double dangling_pr;
    double alpha = DEFAULT_ALPHA;
    vector<double> old_pr;
    vector<size_t>::iterator ci;

    sum_pr = 0;
    dangling_pr = 0;

    // string file_vector = "vector" + to_string(my_rank) +".txt";
    // string file_matrix = "mat" + to_string(my_rank) +".txt";
    string file_vector = "vector0.txt";
    string file_matrix = "mat0.txt";

    vector< vector<size_t> > rows;
    vector<size_t> num_outgoing;
    int len;

    read_vector(file_vector,len, num_outgoing);
    read_matrix(file_matrix,len,rows);

    int num_rows = len;
    pr.resize(len*proc_n);
    old_pr.resize(len);
    pr[0] = 1;

    int num_iterations = 0;
    double diff = 1;
    int exit_flag = 0;

    vector<double> temp_pr;
    temp_pr.resize(len);
    /********************************/

    sprintf(filename, "output%d.txt", my_rank);
    f = fopen(filename, "w");

    MPI_Barrier(MPI_COMM_WORLD);

    while (diff > convergence && num_iterations < max_iterations){
        if (my_rank == 0) {
            for (i = 1; i < proc_n; i++)
                MPI_Isend(&pr[0], len*proc_n, MPI_INT, i, tag, MPI_COMM_WORLD, &request);
        } else {
            // t1 = MPI_Wtime();
            MPI_Irecv(&pr[0], len*proc_n, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
            if(pr[0] == -1){
                exit_flag = 1;
            }
            MPI_Wait(&request, &status);
            // t2 = MPI_Wtime();
            // fprintf(f, "[Worker %d] Time to receive from master is %f\n", my_rank, t2-t1);
        }
        // MPI_Barrier(MPI_COMM_WORLD);
        if(exit_flag) break;

        if(my_rank > 0){

        // if(pr[0] == -1) { break; }

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
                pr[i +my_rank*len] = h + one_Av + one_Iv;
                diff += fabs(pr[i +my_rank*len] - old_pr[i +my_rank*len]);
                // printf("%f\n",dangling_pr);
                // if(num_iterations == 12){
                //     printf("h: %f one_Av: %f one_Iv: %f\n",h, one_Av, one_Iv);        
                // }
            }
            // printf("diff: %f\n",diff);
            // if (trace) {
                // cout << num_iterations << ": ";
            // }

        }

        if (my_rank == 0) {
            // t1 = MPI_Wtime();
            for (i = 0; i < proc_n - 1; i++){
                MPI_Irecv(&temp_pr[0], len, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[i]);
                MPI_Wait(&requests[i], &status);
                for(int ii=0; ii< len; ii++){ //
                    pr[(status.MPI_SOURCE)*len + ii] = temp_pr[ii];
                    if(temp_pr[ii] == -1){
                        pr[0] = -1;
                        exit_flag = 1;
                    }
                }
            }
            // t2 = MPI_Wtime();
            // fprintf(f, "[Master] Time to receive from workers is %f\n", t2-t1);
        } else {
            // t1 = MPI_Wtime();
            MPI_Isend(&pr[my_rank*len], len, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            // t2 = MPI_Wtime();
            // fprintf(f, "[Worker %d] Time to send to master is %f\n", my_rank, t2-t1);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        num_iterations++;
    }
    pr[0] = -1;
    if(my_rank > 0){
        MPI_Isend(&pr[0], len, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
    }

    fclose(f);

    MPI_Finalize();
    return 0;
}

