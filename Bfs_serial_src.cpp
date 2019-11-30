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
const double DEFAULT_CONVERGENCE = 0.001;
const unsigned long DEFAULT_MAX_ITERATIONS = 10000;
int trace = 1;

void read_matrix(const string &filename, size_t *edgesDest, size_t *offsets, size_t num_vertices){
    istream *infile;
    infile = new ifstream(filename.c_str());
    string line;

    size_t count = 0;
    size_t prev = -1;
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
    }
    prev++;
    for (; prev <= num_vertices; prev++)
        offsets[prev] = count;
}

int main(){

    struct timespec start, stop; 
    struct timespec start_tot, stop_tot; 
    double time;

    double convergence = DEFAULT_CONVERGENCE;
    unsigned long max_iterations = DEFAULT_MAX_ITERATIONS;

    double sum_pr;
    double dangling_pr;
    double alpha = DEFAULT_ALPHA;
    int i;

    sum_pr = 0;
    dangling_pr = 0;

    string file_metadata = "/staging/vkp2/tye69227/256/graph.metadata";
    string file_matrix = "/staging/vkp2/tye69227/256/graph.txt";

    size_t num_vertices;
    size_t num_edges;

    FILE *fp = fopen(file_metadata.c_str(), "r");
    fscanf(fp, "%lu %lu", &num_vertices, &num_edges);
    fclose(fp);

    size_t *edgesDest, *offsets;
    size_t *visited;
    unsigned *labels;

    visited = (size_t *) calloc(num_vertices, sizeof(size_t));
    labels = (unsigned *) calloc(num_vertices, sizeof(unsigned));

    edgesDest = (size_t *) malloc(num_edges * sizeof(size_t));
    offsets = (size_t *) malloc((num_vertices + 1) * sizeof(size_t));

    read_matrix(file_matrix, edgesDest, offsets, num_vertices);

    int num_iterations = 0;
    size_t length_active, new_length_active;


    for (size_t k = 0; k < 1024; k++) {
        size_t s = num_vertices / 1024 * k;
        labels[s] = 1;
        visited[s] = 1;
    }
    length_active = 1024;

/*
    size_t s = 0;
    labels[s] = 1;
    active_vertices[0] = s;
    visited[s] = 1;
    length_active = 1;
*/
    if( clock_gettime( CLOCK_REALTIME, &start_tot) == -1 ) { perror( "clock gettime" );}

    while (num_iterations < max_iterations) {
        if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}

        if (length_active == 0) {
            if( clock_gettime( CLOCK_REALTIME, &stop_tot) == -1 ) { perror( "clock gettime" );}       
            time = (stop_tot.tv_sec - start_tot.tv_sec)+ (double)(stop_tot.tv_nsec - start_tot.tv_nsec)/1e9;
            printf("Total Execution Time: %f sec\n",time);
            int count = 0;
            for (i = 0; i < num_vertices; i++)
                if (labels[i] == 0) {
                    count++;
                }
            cout << "Unvisited: " << count;
            cout << endl;
            break;                 
        }

        new_length_active = 0;
        for (i = 0; i < num_vertices; i++) {
            if (labels[i] != num_iterations + 1)
                continue;
            size_t l = offsets[i];
            size_t r = offsets[i+1];

            for (size_t ci = l; ci < r; ci++) {
                size_t vid = edgesDest[ci];
                if (visited[vid] == 0) {
                    visited[vid] = 1;
                    //active_vertices[new_length_active] = vid;
                    labels[vid] = num_iterations + 2;
                    new_length_active++;
                }
            }
        }

        length_active = new_length_active;
        num_iterations++;
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}       
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        printf("Iter Time: %f sec\n",time);
    }
}

