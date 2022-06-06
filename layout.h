#ifndef LAYOUT_H
#define LAYOUT_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

typedef struct
{
    uint32_t ver_n, ver_m;   // number of edges in list
    uint32_t *vex;           // vertex list
    uint32_t *read_sta;      // sta pos on read
    uint32_t *read_end;      // end pos on read
    int32_t *edge_sta;      // sta pos on edge seq
    int32_t *edge_end;      // end pos on edge seq
}EdgeList_t;

// bi-directed graph
typedef struct
{
    uint32_t adjvex;        // id of next/pre vertex, id<<1 + 0/1 B/E
    uint32_t reduce;        // 0:exist 1:reduced
    uint32_t sta_pos;       // pos of unmatched seq
    uint32_t end_pos;       // pos of unmatched seq
    uint32_t id;            // id of edges of v->w
    uint32_t str;
    uint32_t utg_id;
    char *edge_seq;         // edge sequence
    EdgeList_t edge_list;   // edge list
} strEdge;

typedef struct
{
    uint32_t arc_n; // number of adjacent next node
    strEdge *out_edge;  // out edge list
    uint32_t arc_m; // number of adjacent pre node
    strEdge *in_edge;   // in edge list
} strVertex;

typedef struct
{
    strVertex *ver; // vertex
    uint32_t nextId; // number of vertices
    uint32_t allocN; // dynamiclly alloc memory size
} StringGraph;

// typedef struct UNITIG
// {
//     uint32_t id;    // id of unitig
//     uint32_t str;
//     uint32_t len;   // length of unitig
//     char *edge_seq; // sequence of unitig
//     uint32_t seq_n;
//     EdgeList_t edge_list;   // edge list
// }UTG_t;

// typedef struct CONTIG
// {
//     uint32_t id;
//     uint32_t str;
//     uint32_t len;
//     uint32_t seqn_n;
//     char *edge_seq;
// }CTG_t;

int layout_workflow(const char *read_fastq, const char *temp_ove_dir, const char *temp_ove_link_dir);

#endif