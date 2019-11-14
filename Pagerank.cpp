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
#define MODE 0 // 0 - ring 1 - s2c2
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

void ring_pagerank(int id, MPI_Status status, int proc_n, int tag){
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
        for(int sas_i = 0; sas_i < proc_n-1; sas_i++){
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
	                old_pr[i] = pr[i]; // old_pr[i] = pr[i] / sum_pr;
	            }
	        }
	        /*
	         * After normalisation the elements of the pagerank vector sum
	         * to one
	         */
	        // sum_pr = 1;
	        
	        /* An element of the A x I vector; all elements are identical */
	        double one_Av = alpha * dangling_pr / num_rows;

	        /* An element of the 1 x I vector; all elements are identical */
	        double one_Iv = (1 - alpha) / num_rows; //  double one_Iv = (1 - alpha) * sum_pr / num_rows;

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
	                h += h_v * (old_pr[*ci]/sum_pr);
	            }
	            h *= alpha;
	            pr[i] = h + one_Av + one_Iv;
	            diff += fabs(pr[i] - old_pr[i]);
	            // printf("%f\n",dangling_pr);
	            // if(num_iterations == 12){
	            //     printf("h: %f one_Av: %f one_Iv: %f\n",h, one_Av, one_Iv);        
	            // }
	        }
	        MPI_Barrier(MPI_COMM_WORLD);
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


/* This is the master */
// int master_io(int id, MPI_Status status, int proc_n, int tag)
// {

    

//     double convergence = DEFAULT_CONVERGENCE;
//     unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

//     double sum_pr;
//     double dangling_pr;
//     double alpha = DEFAULT_ALPHA;
//     vector<size_t>::iterator ci;
//     int i;

//     sum_pr = 0;
//     dangling_pr = 0;

//     string file_vector = "vector.txt";


//     vector<size_t> num_outgoing;
//     int len;
//     read_vector(file_vector,len, num_outgoing);



//     pr[0] = 1;

//     int num_iterations = 0;
//     double diff = 1;
//     int terminate = 0;
//     int temp_teminate = -1;

//     // while(terminate == 0){
//         for(int k=0; k < proc_n; k++){
//             MPI_Send (&pr[0], len*proc_n, MPI_INT, k, tag, MPI_COMM_WORLD);
//         }
//         // if(temp_teminate == -1){
//         //     terminate = 0;
//         //     break;
//         // }

//         for(int k=0; k < proc_n; k++){
//             MPI_Recv (&temp_pr[0], len, MPI_INT, k, tag, MPI_COMM_WORLD, &status);
//             if(temp_pr[0] == -1){
//                 temp_teminate = 0;
//             }
//             for(int ii=0; ii< len; ii++){
//                 pr[k*len + ii] = temp_pr[ii];
//             }
//         }
//     // }



// }

void s2c2_pagerank(int id, MPI_Status status, int proc_n, int tag, MPI_Request request, MPI_Request requests[])
{
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
    pr.resize(len*proc_n);
    old_pr.resize(len);
    pr[0] = 1;

    int num_iterations = 0;
    double diff = 1;

    vector<double> temp_pr;
    temp_pr.resize(len);

    // iter_pagerank(vector< vector<size_t> > rows, vector<size_t> num_outgoing, vector<double> &pr, double &diff, int &num_iterations, int num_rows)
    // while (diff > convergence && num_iterations < max_iterations){
        // printf("id: %d\n",id);

        // printf("id: %d\n",id);
        if(id == 0){


            for(int k=1; k < proc_n; k++){
                MPI_Isend (&pr[0], len*proc_n, MPI_INT, k, tag, MPI_COMM_WORLD, &request);
            }


        }
        else{
            MPI_Irecv (&pr[0], len*proc_n, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
        }
        if(id > 0){

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
                pr[i +id*len] = h + one_Av + one_Iv;
                diff += fabs(pr[i +id*len] - old_pr[i +id*len]);
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

        if(id == 0){
            for(int k=1; k < proc_n; k++){
                MPI_Irecv (&temp_pr[0], len, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,  &requests[i]);

                MPI_Wait(&requests[k], &status);
                for(int ii=0; ii< len; ii++){
                    pr[(k-1)*len + ii] = temp_pr[ii];
                }
            }

        }
        else{
        	 MPI_Isend (&pr[0], len, MPI_INT, 0, tag, MPI_COMM_WORLD,&request); 
        	 MPI_Wait(&request, &status);
        }


    
    // }

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

    // pr[0] = -1;
    // if(id > 0)
    // 	MPI_Send (&pr[0], len, MPI_INT, 0, tag, MPI_COMM_WORLD);


}

#if MODE == 0
int main(){

    int tag = 50;
    MPI_Status status;  

    MPI_Init (NULL , NULL);
    int my_rank, proc_n;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);
    // printf("%d\n",proc_n);

    //ring model
    ring_pagerank(my_rank,status,proc_n,tag);

    MPI_Finalize();


    // vector<size_t>::iterator ci;
    // for(int i=0; i<x; i++){
    //     for (ci = rows[i].begin(); ci != rows[i].end(); ci++)
    //         if(*ci != 0)
    //             printf("num: %f %d\n",*ci, i);
    // }
    return 0;
}

#else

// int main(){


//     int tag = 50;
//     MPI_Status status;  
// 	MPI_Request request;
//     MPI_Request requests[32];

//     MPI_Init (NULL , NULL);
//     int my_rank, proc_n;
//     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &proc_n);
//     // printf("%d\n",proc_n);

//     //master/worker
//     // MPI_Comm_split( MPI_COMM_WORLD, my_rank == 0, 0, &new_comm);

//     s2c2_pagerank(my_rank, status, proc_n, tag, request, requests);


//     MPI_Finalize();


//     // vector<size_t>::iterator ci;
//     // for(int i=0; i<x; i++){
//     //     for (ci = rows[i].begin(); ci != rows[i].end(); ci++)
//     //         if(*ci != 0)
//     //             printf("num: %f %d\n",*ci, i);
//     // }
//     return 0;

// }

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

    MSIZE = atoi(argv[1]);

    message = (int *) malloc(sizeof(int) * MSIZE);

    message[0] = my_rank * 3;
    message[1] = my_rank * 3 + 1;
    message[2] = my_rank * 3 + 2;

    sprintf(filename, "output%d.txt", my_rank);
    f = fopen(filename, "w");

    MPI_Barrier(MPI_COMM_WORLD);

    if (my_rank == 0)
        printf("[Master] My rank is %d, Total # of proc is %d, Data is %d %d %d\n", my_rank, proc_n, message[0], message[1], message[2]);
    else
        printf("[Worker] My rank is %d, Total # of proc is %d, Data is %d %d %d\n", my_rank, proc_n, message[0], message[1], message[2]);

    for (iteration = 0; iteration < 5; iteration++) {
        if (my_rank == 0) {
            for (i = 1; i < proc_n; i++)
                MPI_Isend(message, MSIZE, MPI_INT, i, tag, MPI_COMM_WORLD, &request);
        } else {
            t1 = MPI_Wtime();
            MPI_Irecv(message, MSIZE, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            t2 = MPI_Wtime();
            printf("[Worker %d] Data is %d %d %d\n", my_rank, message[0], message[1], message[2]);
            fprintf(f, "[Worker %d] Time to receive from master is %f\n", my_rank, t2-t1);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if (my_rank == 0) {
            t1 = MPI_Wtime();
            for (i = 0; i < proc_n - 1; i++)
                MPI_Irecv(&message[MSIZE/p*i], MSIZE/p, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[i]);
            for (i = 0; i < proc_n - 1; i++)
                MPI_Wait(&requests[i], &status);
            t2 = MPI_Wtime();
            fprintf(f, "[Master] Time to receive from workers is %f\n", t2-t1);
        } else {
            t1 = MPI_Wtime();
            MPI_Isend(&message[MSIZE/p*(my_rank - 1)], MSIZE/p, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            t2 = MPI_Wtime();
            fprintf(f, "[Worker %d] Time to send to master is %f\n", my_rank, t2-t1);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    fclose(f);

    MPI_Finalize();
    return 0;
}


#endif


