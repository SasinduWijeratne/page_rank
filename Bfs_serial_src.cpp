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
const unsigned long DEFAULT_MAX_ITERATIONS = 100; //10000;
int trace = 1;

void read_matrix(const string &filename, size_t *edgesDest, size_t *offsets){
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
}

int main(){

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

    size_t *edgesDest, *offsets;
    size_t *active_vertices, *visited;
    unsigned *labels;

    active_vertices = (size_t *) malloc(num_vertices * sizeof(size_t));
    visited = (size_t *) calloc(num_vertices, sizeof(size_t));
    labels = (unsigned *) calloc(num_vertices, sizeof(unsigned));

    edgesDest = (size_t *) malloc(num_edges * sizeof(size_t));
    offsets = (size_t *) malloc((num_vertices + 1) * sizeof(size_t));

    offsets[num_vertices] = num_edges;
    read_matrix(file_matrix, edgesDest, offsets);
    
    int num_iterations = 0;
    size_t length_active, new_length_active;

    for (size_t k = 0; k < 1024; k++) {
        size_t s = num_vertices / 1024 * k;
        labels[s] = 1;
        active_vertices[k] = s;
        visited[s] = 1;
    }
    length_active = 1024;

    if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
    // iter_pagerank(vector< vector<size_t> > rows, vector<size_t> num_outgoing, vector<double> &pr, double &diff, int &num_iterations, int num_rows)
    while (num_iterations < max_iterations) {

        if (length_active == 0) {
            int count = 0;
            for (i = 0; i < num_vertices; i++)
                if (labels[i] == 0)
                    count++;
            cout << "Unvisited: " << count;
            cout << endl;
            break;                 
        }

        new_length_active = 0;
        for (i = 0; i < length_active; i++) {
            for (size_t ci = offsets[i]; ci < offsets[i+1]; ci++) {
                size_t vid = edgesDest[ci];
                if (!visited[vid]) {
                    visited[vid] = 1;
                    active_vertices[new_length_active] = vid;
                    labels[vid] = num_iterations + 1;
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

