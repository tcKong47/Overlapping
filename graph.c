#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "main.h"
#include "ktime.h"
#include "overlapping.h"
#include "bit_operation.h"

// #define _DEBUG
// #define _DDEBUG

void init_graph(Graph *Dp_graph, uint32_t max_v_num)
{
	// Dp_graph = (Graph* )calloc(1,sizeof(Graph));
	// if(Dp_graph == NULL) fprintf(stderr, "[Wrong] Failed to allocate graph memory!\n");
	Dp_graph->vnode = (VNode *)calloc(max_v_num, sizeof(VNode));
	if (Dp_graph->vnode == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB Dp_graph->vnode memory error!\n", __func__, max_v_num * sizeof(VNode) / 1024 / 1024 / 1024);
		exit(1);
	}
}

void free_graph(Graph *Dp_graph)
{
	if (Dp_graph->vnode != NULL){free(Dp_graph->vnode);Dp_graph->vnode = NULL;}
	// if (Dp_graph != NULL)	free(Dp_graph);
}

void find_longest_path(Graph *graph, uint32_t vi, PATH_t *dist_path, uint32_t vertex_num)
{
	float temp = 0;
    float penalty;
	float current_dist = 0;
	int32_t pre_node = -1;
    int weight, dir;
	// int cur_dir = -1;
	uint32_t adjvex;

    //travel the processer of vertex vi
    pANode arcnode;
	int i;
	for (i = 0; i < graph->vnode[vi].anode_num; ++i)
	{
		arcnode = &graph->vnode[vi].preedge[i];
		weight = arcnode->weight;
		penalty = arcnode->penalty;
		dir = arcnode->dir;
		adjvex = arcnode->adjvex;
		temp = dist_path[adjvex].dist + weight - penalty;

		// if (dist_path[adjvex].dir == -1) //deal with the start node
		// {
		// 	dist_path[adjvex].dir = dir;
		// }

		if (current_dist <= temp && dist_path[adjvex].dir == dir) //only walk to the same direction
        {
            current_dist = temp;
            pre_node = adjvex;
            // cur_dir = dir;
        }
	}

	// dist_path[vi].dir = cur_dir; //record the current direction
	dist_path[vi].dist = current_dist; // change the distance because there are overlaps between mems
	dist_path[vi].pre_node = pre_node; //the front node

	if (vi >= vertex_num) printf("[%s] use unlegal memory, %d-%d, something wrong with it...", __func__, vi, vertex_num);
	return;
}

int sparse_dynamic_programming(Graph *graph, uint32_t vertex_num, MR_t *vertex_mr, PATH_t *dist_path, uint8_t *out_degree)
{
	uint32_t i;

	for (i = 0; i < vertex_num; ++i)
	{
		if (graph->vnode[i].anode_num > 0)
		{
			find_longest_path(graph, i, dist_path, vertex_num);
		}
	}

	#ifdef _DEBUG
	int j;
	printf("show the path");
	for (i = 0; i < vertex_num; ++i)
	{
		j = i;
		if ((out_degree[i] == 0))
		{
			printf("\nvertex: %d\tdist = %f\n", i, dist_path[i].dist);
			printf("path: ");
			printf("%d(%d,%ld,%ld,%ld,%ld)->", i, vertex_mr[i].t_id, vertex_mr[i].ts, vertex_mr[i].te, vertex_mr[i].qs, vertex_mr[i].qe);

			j = dist_path[j].pre_node;
			while (j != -1)
			{
				printf("%d(%d,%ld,%ld,%ld,%ld)->", j, vertex_mr[j].t_id, vertex_mr[j].ts, vertex_mr[j].te, vertex_mr[j].qs, vertex_mr[j].qe);
				j = dist_path[j].pre_node;
			}
		}
	}
	printf("\nend\n");
	#endif

	return 0;
}

int create_graph(thread_buf_t *buf)
{
	int32_t i, j, weight, gap;
	int32_t ove_q, ove_t;
	uint32_t non_ioslated_point = 0;
	uint32_t thre_num;
	uint32_t vertex_num = buf->v_n;
	int search_step = (vertex_num < 10)? vertex_num : 10, t_id, index_n;
	uint32_t *max_dis = NULL;
	uint8_t *out_degree = NULL;
	Graph *graph = &buf->Dp_graph;
	MR_t *vertex_mr = buf->vertex_mr;
	PATH_t *dist_path = buf->path;
	uint32_t *max_index = buf->max_index;

#ifdef _DDEBUG
	printf("qid\thid\tcov\tqstrand\tqs\tqe\ttid\ttstrand\tts\tte\tlen\n");
	for (int vi = 0; vi < vertex_num; ++vi)
	{
		printf("%d\t%d\t%d\t%ld\t%ld\t%d\t%d\t%ld\t%ld\n", vi, vertex_mr[vi].cov, vertex_mr[vi].qstr, vertex_mr[vi].qs, vertex_mr[vi].qe, vertex_mr[vi].t_id, vertex_mr[vi].tstr, vertex_mr[vi].ts, vertex_mr[vi].te);
	}
	printf("\n");
#endif

	// re inital
	for (i = 0; i < vertex_num; ++i)
	{
		graph->vnode[i].anode_num = 0;
		graph->vnode[i].preedge = (ANode*)calloc(search_step, sizeof(ANode));
		if (graph->vnode[i].preedge == NULL)
		{
			fprintf(stderr, "[%s] calloc %ldGB graph->vnode[i].preedge memory error!\n", __func__, search_step * sizeof(ANode) / 1024 / 1024 / 1024);
			exit(1);
		}
		dist_path[i].dir = vertex_mr[i].qstr == vertex_mr[i].tstr ? 0 : 1;
		dist_path[i].dist = vertex_mr[i].cov;
		dist_path[i].pre_node = -1;
	}

	max_dis = (uint32_t *)calloc(vertex_num, sizeof(uint32_t));
	if (max_dis == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB max_dis memory error!\n", __func__, vertex_num * sizeof(uint32_t) / 1024 / 1024 / 1024);
		exit(1);
	}

	out_degree = (uint8_t *)calloc(vertex_num, sizeof(uint8_t));
	if (out_degree == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB out_degree memory error!\n", __func__, vertex_num * sizeof(uint8_t) / 1024 / 1024 / 1024);
		exit(1);
	}

	for (i = 0; i < vertex_num - 1; ++i)
	{
		thre_num = (vertex_num < (i + search_step))? vertex_num : (i + search_step);
		for (j = i + 1; j < thre_num; ++j)
		{
			if (vertex_mr[i].t_id != vertex_mr[j].t_id)	break;
			if ((vertex_mr[i].tstr ^ vertex_mr[j].tstr) != (vertex_mr[i].qstr ^ vertex_mr[j].qstr))	continue;

			ove_q = vertex_mr[j].qs - vertex_mr[i].qe - 1;
			ove_t = vertex_mr[j].tstr == vertex_mr[j].qstr ? (int32_t)(vertex_mr[j].ts - vertex_mr[i].te - 1) : (int32_t)(vertex_mr[i].ts - vertex_mr[j].te - 1);

			gap = (int32_t)abs(ove_q - ove_t);
			// printf("%d->%d, %d %d %d\n",i,j,ove_q,ove_t,gap);

			// if (ove_q > -k && ove_t > -k && gap < fmin(abs(ove_t), abs(ove_q)))  //5  1/error rate
			if (ove_q > -k && ove_t > -k && (gap <= abs(ove_q)/5 || gap < k))
            {
            	graph->vnode[j].preedge[graph->vnode[j].anode_num].adjvex = i;
				weight = ove_q >= 0 ? vertex_mr[j].cov: vertex_mr[j].cov + ove_q;
				graph->vnode[j].preedge[graph->vnode[j].anode_num].weight = weight;
				graph->vnode[j].preedge[graph->vnode[j].anode_num].penalty = gap / (float)weight;
				graph->vnode[j].preedge[graph->vnode[j].anode_num++].dir = vertex_mr[j].tstr == vertex_mr[j].qstr ? 0 : 1;
				out_degree[i] = 1;

        		non_ioslated_point++;
				// printf("%d->%d, weight %d, penalty %.3f\n",i,j,weight,graph->vnode[j].preedge[graph->vnode[j].anode_num-1].penalty);
            }
			if (j >= vertex_num) printf("[%s] 1 use unlegal memory, %d-%d, something wrong with it...", __func__, j, vertex_num);
			if (graph->vnode[j].anode_num > search_step) printf("[%s] use unlegal memory, %d-%d, something wrong with it...", __func__, graph->vnode[j].anode_num, search_step);
		}
	}

	// if (non_ioslated_point != 0)
	{
		sparse_dynamic_programming(graph, vertex_num, vertex_mr, dist_path, out_degree);
	}

	index_n = 0;
	t_id = vertex_mr[0].t_id;
	for (i = 0; i < vertex_num; ++i)
	{
		if (vertex_mr[i].t_id != t_id)
		{
			if (dist_path[max_index[index_n]].pre_node != -1 || max_dis[index_n] > k)	index_n++;
			t_id = vertex_mr[i].t_id;
			max_index[index_n] = i;
		}
		if (max_dis[index_n] < dist_path[i].dist)
		{
			max_dis[index_n] = dist_path[i].dist;
			max_index[index_n] = i;
		}
	}

	if (dist_path[max_index[index_n]].pre_node != -1 || max_dis[index_n] > k) //deal with the last node
		index_n++;

	buf->max_index_n = index_n;
	if (index_n > vertex_num) printf("[%s] 2 use unlegal memory, %d-%d, something wrong with it...", __func__, index_n, vertex_num);

#ifdef _DEBUG
	for (size_t i = 0; i < index_n; i++)
	{
		printf("%d %d %d\n", vertex_mr[i].t_id, max_index[i], max_dis[i]);
	}
#endif

	for (i = 0; i < vertex_num; ++i)
	{
		if (graph->vnode[i].preedge != NULL) {free(graph->vnode[i].preedge);graph->vnode[i].preedge = NULL;}
	}
	if (max_dis != NULL) {free(max_dis);max_dis = NULL;}
	if (out_degree != NULL) {free(out_degree);out_degree = NULL;}

	return 0;
}