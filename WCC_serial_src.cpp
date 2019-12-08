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

/*mpic++ your_code_file.c
Execution

mpirun -np <no. of Processors> ./a.out
*/

const double DEFAULT_ALPHA = 0.85;
const double DEFAULT_CONVERGENCE = 0.001;
const unsigned long DEFAULT_MAX_ITERATIONS = 100; //10000;
int trace = 1;

void read_matrix(const string &filename, size_t num_vertices, size_t num_edges, size_t *edgesDest, size_t *offsets, size_t *edgesDest2, size_t *offsets2){
    istream *infile;
    infile = new ifstream(filename.c_str());
    string line;

    size_t count = 0;
    size_t count2 = 0;
    size_t prev = -1;
    size_t prev2 = -1;

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

        if (prev2 != dest) {
            prev2++;
            for (; prev2 < dest; prev2++)
                offsets2[prev2] = count2;
            offsets2[dest] = count2;
            prev2 = dest;
        }
        edgesDest2[count2] = src;
        count2++;
    }

    prev++;
    for (; prev <= num_vertices; prev++)
        offsets[prev] = num_edges;

    prev2++;
    for (; prev2 <= num_vertices; prev2++)
        offsets2[prev2] = num_edges;

}

void wcc(){

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

    double *pr, *rank;
    size_t *edgesDest, *offsets;
    size_t *edgesDest2, *offsets2;

    pr = (double *) malloc(num_vertices*sizeof(size_t));
    rank = (double *) malloc(num_vertices*sizeof(size_t));

    edgesDest = (size_t *) malloc(num_edges * sizeof(size_t));
    offsets = (size_t *) malloc((num_vertices + 1) * sizeof(size_t));

    edgesDest2 = (size_t *) malloc(num_edges * sizeof(size_t));
    offsets2 = (size_t *) malloc((num_vertices + 1) * sizeof(size_t));

    read_matrix(file_matrix, num_vertices, num_edges, edgesDest, offsets, edgesDest2, offsets2);
    
    for (size_t k = 0; k < num_vertices; k++)
        pr[k] = num_vertices; // initialize pr[:] = INF

    for (size_t k = 0; k < num_vertices; k++)
        rank[k] = k;

    int num_iterations = 0;
    int hasUpdate = 1;

    if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
    while (num_iterations < max_iterations && hasUpdate) {
        hasUpdate = 0;
        for (size_t v = 0; v < num_vertices; v++) {
            for (size_t ci = offsets2[v]; ci < offsets2[v+1]; ci++) {
                size_t u = edgesDest2[ci];
                if (pr[u] < rank[v]) {
                    rank[v] = pr[u];
                    hasUpdate = 1;
                }
            }
            pr[v] = rank[v];
        }
        num_iterations++;
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}       
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("Round 1 Iter %d Time: %f sec\n", num_iterations, time);
    }

    num_iterations = 0;
    hasUpdate = 1;
    if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
    while (num_iterations < max_iterations && hasUpdate) {
        hasUpdate = 0;
        for (size_t v = 0; v < num_vertices; v++) {
            for (size_t ci = offsets[v]; ci < offsets[v+1]; ci++) {
                size_t u = edgesDest[ci];
                if (pr[u] < rank[v]) {
                    rank[v] = pr[u];
                    hasUpdate = 1;
                }
            }
            pr[v] = rank[v];
        }
        num_iterations++;
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}       
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("Round 2 Iter %d Time: %f sec\n", num_iterations, time);
    }   

}

int main(){
    wcc();
    return 0;
}