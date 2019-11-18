#include <stdlib.h>
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
using namespace std;
#define MODE 0 // 0 - ring 1 - s2c2
/*mpic++ your_code_file.c
Execution

mpirun -np <no. of Processors> ./a.out
*/

const double DEFAULT_ALPHA = 0.85;
const double DEFAULT_CONVERGENCE = 0.00001;
const unsigned long DEFAULT_MAX_ITERATIONS = 100; //10000;
int trace = 1;

void read_matrix(const string &filename, size_t *edgesDest, size_t *offsets, size_t *recip_offsets, int num_vertices){
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

void pagerank(){

    struct timespec start, stop; 
    double time;

    double convergence = DEFAULT_CONVERGENCE;
    unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

    double sum_pr;
    double dangling_pr;
    double alpha = DEFAULT_ALPHA;
    int i;

    sum_pr = 0;
    dangling_pr = 0;

    string file_metadata = "graph.metadata";
    string file_matrix = "graph.txt";

    size_t num_vertices;
    size_t num_edges;

    FILE *fp = fopen(file_metadata.c_str(), "r");
    fscanf(fp, "%lu %lu", &num_vertices, &num_edges);
    fclose(fp);

    double *pr, *old_pr;
    size_t *edgesDest, *offsets, *recip_offsets;
    pr = (double *) calloc(num_vertices , sizeof(double));
    old_pr = (double *) calloc(num_vertices , sizeof(double));
    edgesDest = (size_t *) malloc(num_edges * sizeof(size_t));
    offsets = (size_t *) malloc((num_vertices + 1) * sizeof(size_t));
    recip_offsets = (size_t *) malloc((num_vertices + 1) * sizeof(size_t));

    read_matrix(file_matrix, edgesDest, offsets, recip_offsets,num_vertices);

    // for(int hh = 0; hh < num_vertices; hh++){
    //     printf("%zu %zu\n",recip_offsets[hh],offsets[hh]);
    // }
    
    pr[0] = 1;

    int num_iterations = 0;
    double diff = 1;
    double pre_diff = -1;

    if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
    // iter_pagerank(vector< vector<size_t> > rows, vector<size_t> num_outgoing, vector<double> &pr, double &diff, int &num_iterations, int num_rows)
    while (diff > convergence && num_iterations < max_iterations && fabs(pre_diff-diff) > 0.0000001 ) {

        sum_pr = 0;
        dangling_pr = 0;
        
        for (size_t k = 0; k < num_vertices; k++) {
            double cpr = pr[k];
            sum_pr += cpr;
            if (offsets[k+1] - offsets[k] == 0) {
                dangling_pr += cpr;
            }
        }

        if (num_iterations == 0) {
            old_pr[0] = 1;
            diff = 999;
        } else {
            /* Normalize so that we start with sum equal to one */
            pre_diff = diff;
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
        double one_Av = alpha * dangling_pr / num_vertices;

        /* An element of the 1 x I vector; all elements are identical */
        double one_Iv = (1 - alpha) * sum_pr / num_vertices;

        /* The difference to be checked for convergence */

        for (i = 0; i < num_vertices; i++) {
            /* The corresponding element of the H multiplication */
            double h = 0.0;
            for (size_t ci = offsets[i]; ci < offsets[i+1]; ci++) {
                /* The current element of the H vector */
                double h_v = (recip_offsets[edgesDest[ci]+1]-recip_offsets[edgesDest[ci]])
                    ? 1.0 / (recip_offsets[edgesDest[ci]+1]-recip_offsets[edgesDest[ci]])
                    : 0.0;
                h += h_v * old_pr[edgesDest[ci]];
            }
            h *= alpha;
            pr[i] = h + one_Av + one_Iv;
//            diff += fabs(pr[i] - old_pr[i]);
        }

/*
        double pr_sum = 0, old_sum = 0; 
        cout << "    pr: ";
        for (i = 0; i < num_vertices; i++) {
            cout << pr[i] << " ";
            pr_sum += pr[i];
        }
        //cout << endl << pr_sum << endl;
        cout << endl << "old_pr: ";
        for (i = 0; i < num_vertices; i++) {
            cout << old_pr[i] << " ";
            old_sum += old_pr[i];
        }
        cout << endl;
        //cout << endl << old_sum << endl;
*/
        num_iterations++;
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}       
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("Iter Time: %f sec\n",time);
        printf("Iter: %d Diff: %f \n\n", num_iterations, diff);
    }


}

#if MODE == 0
int main(){

    // printf("%d\n",proc_n);

    //ring model
    pagerank();


    return 0;
}

#else

// int main(){


//     int tag = 50;
//     MPI_Status status;  
//  MPI_Request request;
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


