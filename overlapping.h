#ifndef OVERLAPPING_H
#define OVERLAPPING_H

#include "main.h"

typedef struct ArcNode
{
	uint32_t adjvex; //id of precursor vertex
	uint32_t weight; //weight of edge
	float penalty;	 //penalty of edge
	int dir;		 //direction of edge
} ANode, *pANode;

typedef struct VertexNode
{
	uint32_t anode_num; //number of adjacent node
	ANode *preedge;		//edge link between current node and precursor
} VNode;

typedef struct dpGraph
{
	VNode *vnode; //list of vertex
} Graph;

typedef struct
{
	uint64_t ts, te;
	uint64_t qs, qe;
	uint32_t t_id;
	int32_t cov;
	uint8_t tstr, qstr;
}MR_t;

typedef struct PATH
{
	int dir;
	float dist;
	int32_t pre_node;
}PATH_t;

typedef struct
{
	uint32_t qid, tid;
	uint32_t ql, tl;
	uint32_t qs, ts;
	uint32_t qe, te;
	uint8_t rev;
	uint32_t mbp;
	uint32_t mln;
	float score;
}OVE_t;

typedef struct
{
	uint64_t n, m, nn; //num of overlapping in collection; simple implematation of symm ove_cl for seed read
	OVE_t *ove; //point to a overlapping collection
}OVE_C_t;

typedef struct
{
	mini_t *seq_mi;
	MR_t *vertex_mr;
	uint32_t v_n, v_m;
	uint32_t max_index_n;
	PATH_t *path;
	uint32_t *max_index;
	Graph Dp_graph;
	uint32_t max_v_num;
	OVE_C_t *ove_cl;
} thread_buf_t;

typedef struct
{
	uint32_t rid; // rid of last read
	uint32_t sta_pos, end_pos; // sta pos in seed read
	uint32_t ts, te; // sta pos in read[rid]
	uint8_t rev;
} cov_data_t;

typedef struct
{
    size_t n, m;
    cov_data_t *a;
} rStack;

typedef struct
{
	uint32_t rid;
	uint32_t *cov;
	uint32_t n;
	double ave_cov;
} read_cov_t;

typedef struct
{
	uint32_t *rlen;
	uint32_t *is_idx;  // 0:indexed read, 1:not indexed read
	uint32_t *not_idx; // 0:read can be indexed, 1:read cannot be indexed
	uint32_t *iter_idx;
	uint32_t iter_idx_n;
	uint32_t *iter_map;
	uint32_t iter_map_n;
}read_stat_t;


uint32_t finding_overlapping(const char *read_fastq, const char *index_fastq, const char *temp_ove_dir, const char *temp_ove_link_dir);

#endif