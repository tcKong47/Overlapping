#ifndef ASSEMBLYGRAPH_H
#define ASSEMBLYGRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "layout.h"

typedef struct
{
    size_t n, m;
    uint32_t *a;
} asgStack;

typedef struct cluster
{
    uint32_t *v;
    uint32_t n;
}Cluster_t;

typedef struct edge_list
{
    uint32_t *v;
    uint32_t *w;
    uint32_t n;
}Edge_List;

void init_assembly_graph(StringGraph *string_graph, uint32_t vertex_num);
void free_assembly_graph(StringGraph *string_graph);
int create_assembly_graph(StringGraph *ol_string_graph, StringGraph *bk_string_graph, OVE_C_t *ove_cl, READ_t *target, UTG_t *utg, const char *out_file_dir);

#endif