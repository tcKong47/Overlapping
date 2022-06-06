#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "overlapping.h"

void init_graph(Graph *graph, uint32_t max_v_num);
void free_graph(Graph *graph);
int create_graph(thread_buf_t *buf);

// int create_graph(step_query_t *step, int64_t i, uint32_t pid);

#endif