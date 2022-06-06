#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <assert.h>

#include "bit_operation.h"
#include "overlapping.h"
#include "assemblygraph.h"
#include "kvec.h"
#include "ksw2.h"

#define _DRAW
#define MAX_UINT32 0xffffffff
#define MAX_SEQ_LEN_KSW 2000000 //TODO: find approriate value
#define EDGE_TR 15

uint32_t utg_id = 1;
//global variant
// StringGraph *string_graph = NULL;

void init_assembly_graph(StringGraph *string_graph, uint32_t vertex_num)
{
    string_graph->nextId = 2 * vertex_num - 1;
    string_graph->allocN = 2 * vertex_num;
    string_graph->ver = (strVertex *)calloc(string_graph->allocN, sizeof(strVertex));
    if (string_graph->ver == NULL)
        printf("[%s]Failed to allocate graph->ver memory!\n", __func__);
    for (int i = 0; i < string_graph->allocN; ++i)
    {
        string_graph->ver[i].arc_n = 0;
        string_graph->ver[i].out_edge = (strEdge *)calloc(EDGE_TR * 2, sizeof(strEdge));
        if (string_graph->ver[i].out_edge == NULL)
            printf("[%s] Failed to allocate graph->ver.out_edge memory!\n", __func__);

        string_graph->ver[i].arc_m = 0;
        string_graph->ver[i].in_edge = (strEdge *)calloc(EDGE_TR * 2, sizeof(strEdge));
        if (string_graph->ver[i].in_edge == NULL)
            printf("[%s] Failed to allocate graph->ver.in_edge memory!\n", __func__);

        for (int j = 0; j < EDGE_TR * 2; j++)
        {
            string_graph->ver[i].out_edge[j].edge_list.ver_n = 0;
            string_graph->ver[i].out_edge[j].edge_list.ver_m = 0;
        }
    }
}

void free_assembly_graph(StringGraph *string_graph)
{
    int i,j;
    for (i = 0; i < string_graph->allocN; ++i)
    {
        if (string_graph->ver[i].out_edge != NULL) free(string_graph->ver[i].out_edge);
        if (string_graph->ver[i].in_edge != NULL) free(string_graph->ver[i].in_edge);
    }

    for (i = 0; i < string_graph->allocN; ++i)
    {
        for (j = 0; j < string_graph->ver[i].arc_n; j++)
        {
            if (string_graph->ver[i].out_edge[j].edge_seq != NULL)
            {
                free(string_graph->ver[i].out_edge[j].edge_seq);
                string_graph->ver[i].out_edge[j].edge_seq = NULL;
            }
        }
    }
    for (i = 0; i < string_graph->allocN; ++i)
    {
        for (j = 0; j < string_graph->ver[i].arc_n; j++)
        {
            if (string_graph->ver[i].out_edge[j].edge_list.ver_m > 0)
            {
                if (string_graph->ver[i].out_edge[j].edge_list.vex != NULL)   free(string_graph->ver[i].out_edge[j].edge_list.vex);
                if (string_graph->ver[i].out_edge[j].edge_list.read_sta != NULL)   free(string_graph->ver[i].out_edge[j].edge_list.read_sta);
                if (string_graph->ver[i].out_edge[j].edge_list.read_end != NULL)   free(string_graph->ver[i].out_edge[j].edge_list.read_end);
                if (string_graph->ver[i].out_edge[j].edge_list.edge_sta != NULL)   free(string_graph->ver[i].out_edge[j].edge_list.edge_sta);
                if (string_graph->ver[i].out_edge[j].edge_list.edge_end != NULL)   free(string_graph->ver[i].out_edge[j].edge_list.edge_end);
            }
        }
    }
    
    if (string_graph->ver != NULL)    free(string_graph->ver);
}

// get index of edge v.vt -> w.wt
int asg_get_idx_v_w_edge(StringGraph *string_graph, uint32_t v, uint32_t *vi, uint32_t w, uint32_t step)
{
    int arc_n = 0;
    for (size_t i = 0; i < string_graph->ver[v].arc_n; i++)
    {
        if (string_graph->ver[v].out_edge[i].adjvex == w)
        {
            if (step == 1 && string_graph->ver[v].out_edge[i].reduce != 1)
            {
                vi[arc_n++] = i;
            }
            else if (step == 2 && string_graph->ver[v].out_edge[i].reduce == 0)
            {
                vi[arc_n++] = i;
            }
            else if (step == 3)
            {
                vi[arc_n++] = i;
            }
        }
    }
    return arc_n;
}

int asg_get_idx_v_w_edge_with_id(StringGraph *string_graph, uint32_t v, uint32_t *vi, uint32_t w, uint32_t id)
{
    int arc_n = 0;
    for (size_t i = 0; i < string_graph->ver[v].arc_n; i++)
    {
        if (string_graph->ver[v].out_edge[i].adjvex == w)
        {
            if (string_graph->ver[v].out_edge[i].id == id && string_graph->ver[v].out_edge[i].reduce != 1)
                vi[arc_n++] = i;
        }
    }
    return arc_n;
}

// get index of pre edge v.vt -> w.wt
int asg_get_idx_v_w_preedge(StringGraph *string_graph, uint32_t v, uint32_t *vi, uint32_t w, uint32_t step)
{
    int arc_n = 0;
    for (size_t i = 0; i < string_graph->ver[w].arc_m; i++)
    {
        if (string_graph->ver[w].in_edge[i].adjvex == v)
        {
            if (step == 1 && string_graph->ver[w].in_edge[i].reduce != 1)
            {
                vi[arc_n++] = i;
            }
            else if (step == 2 && string_graph->ver[w].in_edge[i].reduce == 0)
            {
                vi[arc_n++] = i;
            }
        }
    }
    return arc_n;
}

int asg_get_idx_v_w_preedge_with_id(StringGraph *string_graph, uint32_t v, uint32_t *vi, uint32_t w, uint32_t id)
{
    int arc_n = 0;
    for (size_t i = 0; i < string_graph->ver[w].arc_m; i++)
    {
        if (string_graph->ver[w].in_edge[i].adjvex == v)
        {
            if (string_graph->ver[w].in_edge[i].reduce != 1 && string_graph->ver[w].in_edge[i].id == id)
                vi[arc_n++] = i;
        }
    }
    return arc_n;
}

int asg_add_edge_v_w(StringGraph *string_graph, uint32_t v, uint32_t w, uint32_t reduce, uint32_t sta_pos, uint32_t end_pos, char *edge_seq, int32_t *in_idx, int32_t *out_idx)
{
    uint32_t *vi, arc_n;
    vi = (uint32_t *)calloc(EDGE_TR * 2, sizeof(uint32_t));
    // realloc
    if (string_graph->nextId == string_graph->allocN)
    {
        printf("nextId %d, allocN %d\n", string_graph->nextId, string_graph->allocN);
        uint32_t sta = string_graph->nextId, edge_alloc = EDGE_TR * 2;
        string_graph->allocN = string_graph->allocN << 1;
        string_graph->ver = (strVertex *)realloc(string_graph->ver ,string_graph->allocN * sizeof(strVertex));
        for (size_t i = sta; i < string_graph->allocN; ++i)
        {
            string_graph->ver[i].arc_n = 0;
            string_graph->ver[i].out_edge = (strEdge *)calloc(edge_alloc, sizeof(strEdge));
            if (string_graph->ver[i].out_edge == NULL)
                printf("[%s] Failed to allocate graph->ver.out_edge memory!\n", __func__);

            string_graph->ver[i].arc_m = 0;
            string_graph->ver[i].in_edge = (strEdge *)calloc(edge_alloc, sizeof(strEdge));
            if (string_graph->ver[i].in_edge == NULL)
                printf("[%s] Failed to allocate graph->ver.in_edge memory!\n", __func__);

            for (size_t j = 0; j < edge_alloc; j++)
            {
                string_graph->ver[i].out_edge[j].edge_list.ver_n = 0;
                string_graph->ver[i].out_edge[j].edge_list.ver_m = 0;
            }
        }
    }

    // already have edge v->w
    arc_n = asg_get_idx_v_w_edge(string_graph, v, vi, w, 3);

    // add out edge v->w for v
    (*out_idx) = string_graph->ver[v].arc_n;
    string_graph->ver[v].out_edge[string_graph->ver[v].arc_n].adjvex = w;
    string_graph->ver[v].out_edge[string_graph->ver[v].arc_n].reduce = reduce;
    string_graph->ver[v].out_edge[string_graph->ver[v].arc_n].sta_pos = sta_pos;
    string_graph->ver[v].out_edge[string_graph->ver[v].arc_n].end_pos = end_pos;
    string_graph->ver[v].out_edge[string_graph->ver[v].arc_n].edge_seq = edge_seq;
    string_graph->ver[v].out_edge[string_graph->ver[v].arc_n].id = arc_n + 1;
    string_graph->ver[v].arc_n++;

    // add in edge v->w for w
    (*in_idx) = string_graph->ver[w].arc_m;
    string_graph->ver[w].in_edge[string_graph->ver[w].arc_m].adjvex = v;
    string_graph->ver[w].in_edge[string_graph->ver[w].arc_m].reduce = reduce;
    string_graph->ver[w].in_edge[string_graph->ver[w].arc_m].sta_pos = sta_pos;
    string_graph->ver[w].in_edge[string_graph->ver[w].arc_m].end_pos = end_pos;
    string_graph->ver[w].in_edge[string_graph->ver[w].arc_m].edge_seq = edge_seq;
    string_graph->ver[w].in_edge[string_graph->ver[w].arc_m].id = arc_n + 1;
    string_graph->ver[w].arc_m++;

    if (vi != 0)    free(vi);
    return 0;
}

int get_v_arc_n(StringGraph *string_graph, uint32_t v, uint32_t *vi)
{
    int arc_n = 0;
    for (size_t i = 0; i < string_graph->ver[v].arc_n; i++)
    {
        if (string_graph->ver[v].out_edge[i].reduce != 1) // edge not reduced
        {
            vi[arc_n++] = i;
        }
    }
    return arc_n;
}

int get_v_arc_m(StringGraph *string_graph, uint32_t v, uint32_t *vi)
{
    int arc_m = 0;
    for (size_t i = 0; i < string_graph->ver[v].arc_m; i++)
    {
        if (string_graph->ver[v].in_edge[i].reduce != 1) // edge not reduced
        {
            vi[arc_m++] = i;
        }
    }
    return arc_m;
}

// v->w, only for one edge between v->w
int reduce_edge_v_w(StringGraph *string_graph, uint32_t v, uint32_t w)
{
    int arc_n, arc_m, idx, idx_pre;
    int edge_alloc = EDGE_TR * 2;
    uint32_t *vi, *wi;
    vi = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
    wi = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));

    arc_n = asg_get_idx_v_w_edge(string_graph, v, vi, w, 1);
    if (arc_n == 1)
    {
        idx = vi[0];
        if (string_graph->ver[v].out_edge[idx].reduce != 1)
        {
            string_graph->ver[v].out_edge[idx].reduce = 1;
            arc_m = asg_get_idx_v_w_preedge(string_graph, v, wi, w, 1);
            assert(arc_m == 1);
            idx_pre = wi[0];
            string_graph->ver[w].in_edge[idx_pre].reduce = 1;
            return idx;
        }
    }
    if (vi != NULL) free(vi);
    if (wi != NULL) free(wi);
    return -1;
}

// v->w, only for one edge between v->w
int reduce_edge_v_w_with_id(StringGraph *string_graph, uint32_t v, uint32_t w, uint32_t id)
{
    int arc_n, arc_m;
    int edge_alloc = EDGE_TR * 2;
    uint32_t *vi, *wi;
    vi = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
    wi = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));

    arc_n = asg_get_idx_v_w_edge(string_graph, v, vi, w, 1);
    for (size_t i = 0; i < arc_n; i++)
    {
        if (string_graph->ver[v].out_edge[vi[i]].reduce != 1 && string_graph->ver[v].out_edge[vi[i]].id == id)
        {
            string_graph->ver[v].out_edge[vi[i]].reduce = 1;
            arc_m = asg_get_idx_v_w_preedge(string_graph, v, wi, w, 1);
            for (size_t j = 0; j < arc_m; j++)
            {
                if (string_graph->ver[w].in_edge[wi[j]].reduce != 1 && string_graph->ver[w].in_edge[wi[j]].id == id)
                {
                    string_graph->ver[w].in_edge[wi[j]].reduce = 1;
                    break;
                }
            }
            break;
        }
    }
    
    if (vi != NULL) free(vi);
    if (wi != NULL) free(wi);
    return -1;
}

int push_edge_list(strEdge *edge, uint32_t vex, uint32_t read_sta, uint32_t read_end, uint32_t edge_sta, uint32_t edge_end)
{
    if (edge->edge_list.ver_n == edge->edge_list.ver_m)
    {
        edge->edge_list.ver_m = edge->edge_list.ver_m == 0 ? 2 : edge->edge_list.ver_m << 1;
        edge->edge_list.vex = (uint32_t *)realloc(edge->edge_list.vex, edge->edge_list.ver_m * sizeof(uint32_t));
        edge->edge_list.read_sta = (uint32_t *)realloc(edge->edge_list.read_sta, edge->edge_list.ver_m * sizeof(uint32_t));
        edge->edge_list.read_end = (uint32_t *)realloc(edge->edge_list.read_end, edge->edge_list.ver_m * sizeof(uint32_t));
        edge->edge_list.edge_sta = (int32_t *)realloc(edge->edge_list.edge_sta, edge->edge_list.ver_m * sizeof(int32_t));
        edge->edge_list.edge_end = (int32_t *)realloc(edge->edge_list.edge_end, edge->edge_list.ver_m * sizeof(int32_t));
    }
    edge->edge_list.vex[edge->edge_list.ver_n] = vex;
    edge->edge_list.read_sta[edge->edge_list.ver_n] = read_sta;
    edge->edge_list.read_end[edge->edge_list.ver_n] = read_end;
    edge->edge_list.edge_sta[edge->edge_list.ver_n] = edge_sta;
    edge->edge_list.edge_end[edge->edge_list.ver_n] = edge_end;
    edge->edge_list.ver_n++;
    return 0;
}

int push_edge_list_at_beg(strEdge *edge, uint32_t vex, uint32_t read_sta, uint32_t read_end, uint32_t edge_sta, uint32_t edge_end)
{
    if (edge->edge_list.ver_n == edge->edge_list.ver_m)
    {
        edge->edge_list.ver_m = edge->edge_list.ver_m == 0 ? 2 : edge->edge_list.ver_m << 1;
        edge->edge_list.vex = (uint32_t *)realloc(edge->edge_list.vex, edge->edge_list.ver_m * sizeof(uint32_t));
        edge->edge_list.read_sta = (uint32_t *)realloc(edge->edge_list.read_sta, edge->edge_list.ver_m * sizeof(uint32_t));
        edge->edge_list.read_end = (uint32_t *)realloc(edge->edge_list.read_end, edge->edge_list.ver_m * sizeof(uint32_t));
        edge->edge_list.edge_sta = (int32_t *)realloc(edge->edge_list.edge_sta, edge->edge_list.ver_m * sizeof(int32_t));
        edge->edge_list.edge_end = (int32_t *)realloc(edge->edge_list.edge_end, edge->edge_list.ver_m * sizeof(int32_t));
    }
    for (int i = edge->edge_list.ver_n-1; i >= 0; i--)
    {
        edge->edge_list.vex[i+1] = edge->edge_list.vex[i];
        edge->edge_list.edge_sta[i+1] = edge->edge_list.edge_sta[i]+edge_end+1;
        edge->edge_list.edge_end[i+1] = edge->edge_list.edge_end[i]+edge_end+1;
        edge->edge_list.read_sta[i+1] = edge->edge_list.read_sta[i];
        edge->edge_list.read_end[i+1] = edge->edge_list.read_end[i];
    }
    edge->edge_list.vex[0] = vex;
    edge->edge_list.read_sta[0] = read_sta;
    edge->edge_list.read_end[0] = read_end;
    edge->edge_list.edge_sta[0] = edge_sta;
    edge->edge_list.edge_end[0] = edge_end;
    edge->edge_list.ver_n++;
    return 0;
}

int olg_transitive_reduction(StringGraph *string_graph, uint32_t step)
{
    uint8_t *mark;
    uint32_t fuzz = 2000; // max gap difference allowed in transitive recuction between two path
    uint32_t longest = 0, n_reduced = 0, edge_tr = EDGE_TR;
    uint32_t vertex_num = string_graph->nextId + 1;
    uint32_t v, w, x, *vi, *wi, arc_n = 0, arc_w = 0;
    mark = (uint8_t *)calloc(vertex_num, sizeof(uint8_t)); //mark[] = 0:vaccant; 1:inplay; 2:eliminate
    vi = (uint32_t *)calloc(2 * edge_tr, sizeof(uint32_t));
    wi = (uint32_t *)calloc(2 * edge_tr, sizeof(uint32_t));

    //sorting edges for each node, checked
    for (v = 0; v < vertex_num; v++)
    {
        if (string_graph->ver[v].arc_n > 1)
            qsort(string_graph->ver[v].out_edge, string_graph->ver[v].arc_n, sizeof(strEdge), compare_arc);
    }

    for (v = 0; v < vertex_num; v++)
    {
        if (step == 0)  arc_n = string_graph->ver[v].arc_n;
        else if (step == 1) arc_n = get_v_arc_n(string_graph, v, vi);
        if (arc_n > 0)
        {
            if (step == 0)
            {
                for (int j = 0; j < arc_n; j++)
                    mark[string_graph->ver[v].out_edge[j].adjvex] = 1;
                longest = (string_graph->ver[v].out_edge[0].end_pos - string_graph->ver[v].out_edge[0].sta_pos + 1) + fuzz;

                for (int j = 0; j < arc_n; j++)
                {
                    w = string_graph->ver[v].out_edge[j].adjvex;
                    if (mark[w] != 1)   continue;
                    arc_w = string_graph->ver[w].arc_n;
                    for (int k = 0; k < arc_w; k++)
                    {
                        x = string_graph->ver[w].out_edge[k].adjvex;
                        if ((string_graph->ver[v].out_edge[j].end_pos - string_graph->ver[v].out_edge[j].sta_pos) + string_graph->ver[w].out_edge[k].end_pos - string_graph->ver[w].out_edge[k].sta_pos <= longest && mark[x] == 1)
                            mark[x] = 2;
                    }
                }
            }
            else if (step == 1)
            {
                for (int j = 0; j < arc_n; j++)
                    mark[string_graph->ver[v].out_edge[vi[j]].adjvex] = 1;
                longest = (string_graph->ver[v].out_edge[vi[0]].end_pos - string_graph->ver[v].out_edge[vi[0]].sta_pos + 1) + fuzz;

                for (int j = 0; j < arc_n; j++)
                {
                    w = string_graph->ver[v].out_edge[vi[j]].adjvex;
                    if (mark[w] != 1)   continue;
                    arc_w = get_v_arc_n(string_graph, w, wi);
                    for (int k = 0; k < arc_w; k++)
                    {
                        x = string_graph->ver[w].out_edge[wi[k]].adjvex;
                        if ((string_graph->ver[v].out_edge[vi[j]].end_pos - string_graph->ver[v].out_edge[vi[j]].sta_pos) + string_graph->ver[w].out_edge[wi[k]].end_pos - string_graph->ver[w].out_edge[wi[k]].sta_pos <= longest && mark[x] == 1)
                            mark[x] = 2;
                    }
                }
            }

            for (int j = 0; j < string_graph->ver[v].arc_n; j++)
            {
                w = string_graph->ver[v].out_edge[j].adjvex;
                if (mark[w] == 2)
                {
                    // reduce v edge v.B->w.B/w.E, and w.B/w.E->v.E
                    int ret = reduce_edge_v_w(string_graph, v, w);
                    if (ret != -1)  n_reduced++;
                    ret = reduce_edge_v_w(string_graph, w ^ 1, v ^ 1);
                    if (ret != -1)  n_reduced++;
                    // printf("remove edge %d->%d\n", v, w);
                }
                mark[w] = 0;
            }
        }
    }

    if (vi != NULL)   free(vi);
    if (wi != NULL)   free(wi);
    if (mark != NULL)   free(mark);
    return n_reduced;
}

int olg_gen_unbranching_path_seq(StringGraph *string_graph, uint32_t *up, uint32_t up_n, READ_t *target, uint32_t step)
{
    int n_reduced = up_n-2, idx, seq_len = 0, sta_pos = 0, end_pos = 0, seq_n = 0, arc_n;
    uint32_t vex, edge_sta = 0, edge_end, read_sta, read_end;
    uint32_t v, w, *vi, *idx_tmp, ver_n;
    uint32_t edge_tr = EDGE_TR;
    int32_t in_idx, out_idx;
    char *seq;
    int seq_size = MAX_SEQ_LEN_KSW;
    vi = (uint32_t *)calloc(edge_tr * 2, sizeof(uint32_t));
    seq = (char *)calloc(seq_size, sizeof(char));
    idx_tmp = (uint32_t *)calloc(up_n, sizeof(uint32_t));
    
    if (step == 0)
    {
        for (uint32_t i = 0; i < up_n-1; i++)
        {
            v = up[i], w = up[i+1];
            // if (up[0]>>1 == 992)
            //     printf("%d->%d,%d",v>>1,w>>1,seq_n);
            arc_n = asg_get_idx_v_w_edge(string_graph, v, vi, w, 1);
            assert(arc_n == 1);
            idx = vi[0];
            idx_tmp[i] = idx;
            for (size_t j = 0; j < string_graph->ver[v].out_edge[idx].edge_list.ver_n; j++)
            {
                vex = string_graph->ver[v].out_edge[idx].edge_list.vex[j];
                sta_pos = string_graph->ver[v].out_edge[idx].edge_list.read_sta[j];
                end_pos = string_graph->ver[v].out_edge[idx].edge_list.read_end[j];
                seq_len = end_pos - sta_pos;
                if ((vex&1) == 0)
                    seq_n = comp_cat_str(seq, seq_n, target[vex>>1].read_seq, sta_pos, seq_len);
                else if ((vex&1) == 1)
                    seq_n = cat_str(seq, seq_n, target[vex>>1].read_seq, sta_pos, seq_len);
                // if (up[0]>>1 == 992)
                //     printf("->%d\t%d[%d-%d]\t",seq_n, vex>>1,sta_pos,end_pos);
            }
            // if (up[0]>>1 == 992)
            //     printf("\n");
            reduce_edge_v_w(string_graph, v, w);
        }
        // printf("add edge %d->%d, %d\n", up[0]>>1, up[up_n - 1]>>1, seq_n);
        asg_add_edge_v_w(string_graph, up[0], up[up_n - 1], 2, 0, seq_n, seq, &in_idx, &out_idx);
        // add edge list
        string_graph->ver[up[0]].out_edge[out_idx].edge_list.ver_n = 0;
        for (uint32_t i = 0; i < up_n-1; i++)
        {
            v = up[i], w = up[i+1];
            for (size_t j = 0; j < string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.ver_n; j++)
            {
                ver_n = string_graph->ver[up[0]].out_edge[out_idx].edge_list.ver_n;
                vex = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.vex[j];
                read_sta = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.read_sta[j];
                read_end = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.read_end[j];
                edge_sta = ver_n == 0 ? 0 : 1 + string_graph->ver[up[0]].out_edge[out_idx].edge_list.edge_end[ver_n - 1];
                seq_len = read_end - read_sta;
                edge_end = edge_sta + seq_len - 1;
                push_edge_list(&string_graph->ver[up[0]].out_edge[out_idx], vex, read_sta, read_end, edge_sta, edge_end);
                push_edge_list(&string_graph->ver[up[up_n-1]].in_edge[in_idx], vex, read_sta, read_end, edge_sta, edge_end);
            }
        }
    }
    else if (step == 1)
    {
        for (uint32_t i = 0; i < up_n - 1; i++)
        {
            v = up[i], w = up[i+1];
            arc_n = asg_get_idx_v_w_edge(string_graph, v, vi, w, 1);
            // printf("%d, %d->%d, %d\n", arc_n, v >> 1, w >> 1, seq_n);
            assert(arc_n == 1);
            idx = vi[0];
            idx_tmp[i] = idx;

            while (seq_n + string_graph->ver[v].out_edge[idx].end_pos > seq_size)  seq_size = seq_size << 1;
            seq = (char *)realloc(seq, seq_size * sizeof(char));
            seq_n = cat_str(seq, seq_n, string_graph->ver[v].out_edge[idx].edge_seq, 0, string_graph->ver[v].out_edge[idx].end_pos);
            // printf("%d:%d %ld\n",seq_n, string_graph->ver[v].out_edge[idx].end_pos, strlen(string_graph->ver[v].out_edge[idx].edge_seq));
            reduce_edge_v_w(string_graph, v, w);
        }
        // printf("add edge %d %d->%d %d, %d\n", up[0]&1, up[0]>>1, up[up_n-1]&1, up[up_n - 1]>>1, seq_n);
        asg_add_edge_v_w(string_graph, up[0], up[up_n - 1], 2, 0, seq_n, seq, &in_idx, &out_idx);
        // add edge list
        uint32_t acc_len = 0;
        string_graph->ver[up[0]].out_edge[out_idx].edge_list.ver_n = 0;
        for (uint32_t i = 0; i < up_n-1; i++)
        {
            v = up[i], w = up[i+1];
            for (size_t j = 0; j < string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.ver_n; j++)
            {
                vex = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.vex[j];
                read_sta = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.read_sta[j];
                read_end = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.read_end[j];
                edge_sta = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.edge_sta[j] + acc_len;
                edge_end = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.edge_end[j] + acc_len;
                push_edge_list(&string_graph->ver[up[0]].out_edge[out_idx], vex, read_sta, read_end, edge_sta, edge_end);
                push_edge_list(&string_graph->ver[up[up_n-1]].in_edge[in_idx], vex, read_sta, read_end, edge_sta, edge_end);
            }
            acc_len += string_graph->ver[v].out_edge[idx_tmp[i]].end_pos;
        }
    }
    else if (step == 2)
    {
        uint32_t newUtgId = 0, newStr = 0;
        for (uint32_t i = 0; i < up_n - 1; i++)
        {
            v = up[i], w = up[i+1];
            arc_n = asg_get_idx_v_w_edge(string_graph, v, vi, w, 1);
            printf("%d, %d->%d, %d\n", arc_n, v >> 1, w >> 1, seq_n);
            assert(arc_n == 1);
            idx = vi[0];
            idx_tmp[i] = idx;
            if (newUtgId < string_graph->ver[v].out_edge[idx].utg_id)  newUtgId = string_graph->ver[v].out_edge[idx].utg_id, newStr = string_graph->ver[v].out_edge[idx].str;
            while (seq_n + string_graph->ver[v].out_edge[idx].end_pos > seq_size)  seq_size = seq_size << 1;
            seq = (char *)realloc(seq, seq_size * sizeof(char));
            if (i == 0)
            {
                seq_n = cat_str(seq, seq_n, string_graph->ver[v].out_edge[idx].edge_seq, 0, string_graph->ver[v].out_edge[idx].end_pos);
                printf("%d:%d %ld\n", seq_n, 0, strlen(seq));
            }
            else if (i >= 1)
            {
                seq_n = cat_str(seq, seq_n, string_graph->ver[v].out_edge[idx].edge_seq, string_graph->ver[v].out_edge[idx].edge_list.edge_end[0], string_graph->ver[v].out_edge[idx].end_pos - string_graph->ver[v].out_edge[idx].edge_list.edge_end[0]);
                printf("%d:%d %ld, %ld\n", seq_n, string_graph->ver[v].out_edge[idx].edge_list.edge_sta[1], strlen(seq), strlen(string_graph->ver[v].out_edge[idx].edge_seq));
            }
            reduce_edge_v_w(string_graph, v, w);
        }
        printf("add edge %d %d->%d %d, %d %ld\n", up[0]&1, up[0]>>1, up[up_n-1]&1, up[up_n - 1]>>1, seq_n, strlen(seq));
        asg_add_edge_v_w(string_graph, up[0], up[up_n - 1], 2, 0, seq_n, seq, &in_idx, &out_idx);
        string_graph->ver[up[0]].out_edge[out_idx].utg_id = newUtgId;
        string_graph->ver[up[0]].out_edge[out_idx].str = newStr;
        string_graph->ver[up[up_n - 1]].in_edge[in_idx].utg_id = newUtgId;
        string_graph->ver[up[up_n - 1]].in_edge[in_idx].str = newStr;

        // add edge list
        uint32_t acc_len = 0;
        string_graph->ver[up[0]].out_edge[out_idx].edge_list.ver_n = 0;
        for (uint32_t i = 0, j = 0; i < up_n - 1; i++)
        {
            v = up[i], w = up[i + 1];
            if (i == 0) j = 0;
            else if (i >= 1)    j = 1;
            for (; j < string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.ver_n; j++)
            {
                vex = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.vex[j];
                read_sta = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.read_sta[j];
                read_end = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.read_end[j];
                edge_sta = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.edge_sta[j] + acc_len;
                edge_end = string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.edge_end[j] + acc_len;
                if (i == 0)
                {
                    push_edge_list(&string_graph->ver[up[0]].out_edge[out_idx], vex, read_sta, read_end, edge_sta, edge_end);
                    push_edge_list(&string_graph->ver[up[up_n - 1]].in_edge[in_idx], vex, read_sta, read_end, edge_sta, edge_end);
                }
                else if (i >= 1)
                {
                    push_edge_list(&string_graph->ver[up[0]].out_edge[out_idx], vex, read_sta, read_end, edge_sta - string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.edge_end[0], edge_end - string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.edge_end[0]);
                    push_edge_list(&string_graph->ver[up[up_n - 1]].in_edge[in_idx], vex, read_sta, read_end, edge_sta - string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.edge_end[0], edge_end - string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.edge_end[0]);
                }
            }
            if (i == 0) acc_len += string_graph->ver[v].out_edge[idx_tmp[i]].end_pos;
            else if (i >= 1)    acc_len += (string_graph->ver[v].out_edge[idx_tmp[i]].end_pos - string_graph->ver[v].out_edge[idx_tmp[i]].edge_list.edge_end[0]);
        }
    }

    // printf("%s\n",seq);
    if (vi != NULL) free(vi);
    if (idx_tmp != NULL) free(idx_tmp);
    return n_reduced;
}

int olg_merge_unbranching_path(StringGraph *string_graph, READ_t *target, uint32_t step)
{
    int debug_i = 0;
    uint32_t cur_v, sta_v = 0, end_v = 0, idx = 0, i, j, vex;
    uint32_t v, *vi, *vi_pre, *mark, w, idx_pre;
    uint32_t *nl, *up;
    uint32_t n_arc, up_n, arc_m, n_merged = 0;
    uint32_t sta_pos, end_pos, seq_len, seq_n = 0;
    uint32_t vertex_num = string_graph->nextId + 1;

    vi = (uint32_t *)calloc(vertex_num, sizeof(uint32_t));
    vi_pre = (uint32_t *)calloc(vertex_num, sizeof(uint32_t));
    mark = (uint32_t *)calloc(vertex_num, sizeof(uint32_t)); // mark visited node
    nl = (uint32_t *)calloc(vertex_num, sizeof(uint32_t));
    up = (uint32_t *)calloc(vertex_num, sizeof(uint32_t));

    for (v = 0; v < vertex_num; v++)
    {
        if (mark[v] == 1)   continue;
        n_arc = 0, mark[v] = 1, up_n = 0;
        nl[up_n++] = v;
        // printf("0 %d\n", v >> 1);

        // move forward
        cur_v = v;
        if (get_v_arc_n(string_graph, cur_v, vi) == 0)    end_v = cur_v;
        while (get_v_arc_n(string_graph, cur_v, vi) == 1 && get_v_arc_m(string_graph, cur_v, vi_pre) <= 1)
        {
            cur_v = string_graph->ver[cur_v].out_edge[vi[0]].adjvex;
            n_arc++;
            mark[cur_v] = 1, end_v = cur_v;
            nl[up_n++] = cur_v;
            // printf("1 %d\n",cur_v>>1);
            if (end_v == v) {sta_v = end_v; goto up_merge;} // circular unitig
        }

        // move backward
        cur_v = v;
        if (get_v_arc_m(string_graph, cur_v, vi) == 0)    sta_v = cur_v;
        while (get_v_arc_m(string_graph, cur_v, vi_pre) == 1 && get_v_arc_n(string_graph, cur_v, vi) <= 1)
        {
            cur_v = string_graph->ver[cur_v].in_edge[vi_pre[0]].adjvex;
            n_arc++;
            mark[cur_v] = 1, sta_v = cur_v;
            // printf("2 %d\n", cur_v >> 1);
            if (sta_v == end_v) goto up_merge;
            nl[up_n++] = cur_v;
        }

up_merge:
        if (n_arc != 0)
        {
            // printf("%d %d, %d %d\n", sta_v&1, sta_v>>1, end_v&1, end_v>>1);
            // printf("%d %d:\t", debug_i, up_n);
            // for (i = 0; i < up_n; i++)   printf("%d->", nl[i]>>1);
            // printf("\n");

            idx = 0, i = up_n - 1, j = 0;
            if (sta_v == end_v) {up[idx++] = sta_v; up_n++;}
            while(nl[i] != end_v)   up[idx++] = nl[i--];
            while(j <= i)   up[idx++] = nl[j++];
            
            // for (i = 0; i < up_n; i++)   printf("%d->", up[i]>>1);
            // printf("\n");
            debug_i++;
            assert(up_n >= 2);
            if (up_n > 2)
            {
                if (step == 0)
                    n_merged += olg_gen_unbranching_path_seq(string_graph, up, up_n, target, step);
                else if (step == 1)
                    n_merged += olg_gen_unbranching_path_seq(string_graph, up, up_n, NULL, step);
                else if (step == 2)
                    n_merged += olg_gen_unbranching_path_seq(string_graph, up, up_n, NULL, step);
            }
            // if (debug_i == 2)   break;
        }
    }

    // deal with the reduce = 0 edge, associate seq on it
    if (step == 0)
    {
        for (v = 0; v < vertex_num; v++)
        {
            for (j = 0; j < string_graph->ver[v].arc_n; j++)
            {
                w = string_graph->ver[v].out_edge[j].adjvex; // edge v->w
                if (string_graph->ver[v].out_edge[j].reduce == 0)
                {
                    char *seq;
                    seq = (char *)calloc(MAX_SEQ_LEN_KSW, sizeof(char));
                    seq_n = 0;
                    // forward edge
                    arc_m = asg_get_idx_v_w_preedge(string_graph, v, vi, w, 2);
                    // printf("arc_m %d %d.%d->%d.%d\n", arc_m, v>>1, v&1, w>>1, w&1);
                    assert(arc_m == 1);
                    idx_pre = vi[0];

                    for (size_t k = 0; k < string_graph->ver[v].out_edge[j].edge_list.ver_n; k++)
                    {
                        vex = string_graph->ver[v].out_edge[j].edge_list.vex[k];
                        sta_pos = string_graph->ver[v].out_edge[j].edge_list.read_sta[k];
                        end_pos = string_graph->ver[v].out_edge[j].edge_list.read_end[k];
                        seq_len = end_pos - sta_pos;
                        // printf("%d[%d-%d]\t",vex>>1,sta_pos,end_pos);
                        if ((vex&1) == 0)
                            seq_n = comp_cat_str(seq, seq_n, target[vex>>1].read_seq, sta_pos, seq_len);
                        else if ((vex&1) == 1)
                            seq_n = cat_str(seq, seq_n, target[vex>>1].read_seq, sta_pos, seq_len);
                    }

                    string_graph->ver[v].out_edge[j].reduce = 2;
                    string_graph->ver[v].out_edge[j].sta_pos = 0;
                    string_graph->ver[v].out_edge[j].end_pos = seq_n;
                    string_graph->ver[v].out_edge[j].edge_seq = seq;

                    string_graph->ver[w].in_edge[idx_pre].reduce = 2;
                    string_graph->ver[w].in_edge[idx_pre].sta_pos = 0;
                    string_graph->ver[w].in_edge[idx_pre].end_pos = seq_n;
                    string_graph->ver[w].in_edge[idx_pre].edge_seq = seq;
                }
            }
        }
    }

    if (vi != NULL)   free(vi);
    if (vi_pre != NULL)   free(vi_pre);
    if (mark != NULL)   free(mark);
    if (nl != NULL)  free(nl);
    if (up != NULL)  free(up);

    return n_merged;
}

double ksw2_extz_align(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape, int *end_D, int *end_I)
{
    int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
    int tl = strlen(tseq), ql = strlen(qseq);
    // printf("%d %d\n",tl,ql);
    uint8_t *ts, *qs, c[256];
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0;
    c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2;
    c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t *)malloc(tl);
    qs = (uint8_t *)malloc(ql);
    for (i = 0; i < tl; ++i)
        ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i)
        qs[i] = c[(uint8_t)qseq[i]];
    ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
    int idt = 0, nidt = 0;
    for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
    {
        // printf("%d%c", ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
        if ((ez.cigar[i] & 0xf) == 0)
            idt = idt + (ez.cigar[i] >> 4);
        else
        {
            nidt = nidt + (ez.cigar[i] >> 4);
        }
    }
    // putchar('\n');
    // printf("score %d\n", ez.score);
    double identity = 1.0 * idt / (idt + nidt);
    // printf("identity=%d/%d=%f\n", idt, idt + nidt, identity);
    // putchar('\n');
    *end_D = *end_I = 0;
    if ("MID"[ez.cigar[ez.n_cigar-1] & 0xf] == 'D') *end_D = ez.cigar[ez.n_cigar-1] >> 4;
    if ("MID"[ez.cigar[ez.n_cigar-1] & 0xf] == 'I') *end_I = ez.cigar[ez.n_cigar-1] >> 4;
    
    // printf("%d%c...%d%c\n", ez.cigar[0] >> 4, "MID"[ez.cigar[0] & 0xf], ez.cigar[ez.n_cigar-1] >> 4, "MID"[ez.cigar[ez.n_cigar-1] & 0xf]);
    free(ez.cigar);
    free(ts);
    free(qs);
    return identity;
}

//reduce tips supported by one read (one arc tip), TODO: reduce tips supported by max_ext reads
int olg_trimming_tips(StringGraph *string_graph, uint32_t tip_tr)
{
    int ret, n_reduced = 0, end_D, end_I;
    double identity;
    uint32_t v, w, x, *vi, *wi, *ei, wx_id;
    uint32_t edge_tr = EDGE_TR;
    uint32_t arc_n, arc_v, arc_x, idx;
    uint32_t seqn_v, seqn_x, seqn_sub_v, seqn_sub_x;
    uint32_t vertex_num = string_graph->nextId + 1;
    vi = (uint32_t *)calloc(edge_tr, sizeof(uint32_t));
    wi = (uint32_t *)calloc(edge_tr, sizeof(uint32_t));
    ei = (uint32_t *)calloc(edge_tr, sizeof(uint32_t));
    char *vseq, *xseq;
    vseq = (char *)calloc(MAX_SEQ_LEN_KSW, sizeof(char));
    xseq = (char *)calloc(MAX_SEQ_LEN_KSW, sizeof(char));

    // check end tip
    for (v = 0; v < vertex_num; v++) // for each node v
    {
        if (get_v_arc_n(string_graph, v, vi) == 0 && get_v_arc_m(string_graph, v, vi) == 0)   continue;
        if (get_v_arc_n(string_graph, v, vi) != 0)   continue; // v has 0 out edge
        if (get_v_arc_m(string_graph, v, vi) != 1)   continue; // v has 1 in edge, w->v
        w = string_graph->ver[v].in_edge[vi[0]].adjvex;
        if (get_v_arc_n(string_graph, w, vi) <= 1) continue; // w has > 1 out edge

        arc_n = asg_get_idx_v_w_edge(string_graph, w, wi, v, 1);
        assert(arc_n == 1);
        idx = wi[0];
        identity = 0;
        if (string_graph->ver[w].out_edge[idx].end_pos > tip_tr)    continue;
        for (size_t j = 0; j < string_graph->ver[w].arc_n; j++)
        {
            if (string_graph->ver[w].out_edge[j].reduce != 2)   continue;
            x = string_graph->ver[w].out_edge[j].adjvex;
            wx_id = string_graph->ver[w].out_edge[j].id;
            if (x == v) continue;
            seqn_v = string_graph->ver[w].out_edge[idx].end_pos;
            seqn_x = string_graph->ver[w].out_edge[j].end_pos;
            seqn_sub_x = seqn_sub_v = seqn_x > seqn_v ? seqn_v : seqn_x;
            // printf("%d->%d:%d %d\n",w>>1,v>>1,seqn_v,seqn_sub_v);
            // printf("%d->%d:%d %d\n",w>>1,x>>1,seqn_x,seqn_sub_x);
            cat_str(xseq, 0, string_graph->ver[w].out_edge[j].edge_seq, 0, seqn_sub_x);
            cat_str(vseq, 0, string_graph->ver[w].out_edge[idx].edge_seq, 0, seqn_sub_v);
            identity = ksw2_extz_align(vseq, xseq, 1, -2, 2, 1, &end_D, &end_I);
            // printf("idt:%f\n", identity);
            if (identity >= 0.9)
            {
                // merge w->v to w->x's edge list
                // uint32_t vex, read_sta, read_end, edge_sta, edge_end, arc_e;
                // for (size_t k = 0; k < string_graph->ver[w].out_edge[idx].ver_n; k++)
                // {
                //     vex = string_graph->ver[w].out_edge[idx].vex[k];
                //     read_sta = string_graph->ver[w].out_edge[idx].read_sta[k];
                //     read_end = string_graph->ver[w].out_edge[idx].read_end[k];
                //     edge_sta = string_graph->ver[w].out_edge[idx].edge_sta[k];
                //     edge_end = string_graph->ver[w].out_edge[idx].edge_end[k];
                    // push_edge_list(&string_graph->ver[w].out_edge[j], vex, read_sta, read_end, edge_sta,edge_end);
                    // arc_e = asg_get_idx_v_w_preedge_with_id(string_graph, w, ei, string_graph->ver[w].out_edge[j].adjvex, string_graph->ver[w].out_edge[j].id);
                    // assert(arc_e == 1);
                    // push_edge_list(&string_graph->ver[string_graph->ver[w].out_edge[j].adjvex].in_edge[ei[0]], vex, read_sta, read_end, edge_sta, edge_end);
                // }
                ret = reduce_edge_v_w(string_graph, w, v);
                assert(ret != -1);
                n_reduced++;

                // TODO: merge v^1->w^1 to x^1->w^1's edge list
                arc_v = asg_get_idx_v_w_edge(string_graph, v^1, vi, w^1, 1);
                assert(arc_v == 1);
                arc_x = asg_get_idx_v_w_edge_with_id(string_graph, x^1, wi, w^1, wx_id);
                assert(arc_x == 1);
                // uint32_t edge_acc_len = 0, edge_len = 0;
                // for (int k = string_graph->ver[v^1].out_edge[vi[0]].ver_n-1; k >= 0; k--)
                // {
                    // vex = string_graph->ver[v^1].out_edge[vi[0]].vex[k];
                    // read_sta = string_graph->ver[v^1].out_edge[vi[0]].read_sta[k];
                    // read_end = string_graph->ver[v^1].out_edge[vi[0]].read_end[k];
                    // edge_sta = string_graph->ver[v^1].out_edge[vi[0]].edge_sta[k];
                    // edge_end = string_graph->ver[v^1].out_edge[vi[0]].edge_end[k];
                    // edge_len = edge_end - edge_sta;
                    // edge_sta = string_graph->ver[x^1].out_edge[wi[0]].end_pos-edge_len-edge_acc_len;
                    // edge_end = string_graph->ver[x^1].out_edge[wi[0]].end_pos-edge_acc_len;
                    // edge_acc_len += edge_len;
                    // push_edge_list(&string_graph->ver[x^1].out_edge[wi[0]], vex, read_sta, read_end, edge_sta, edge_end);
                    // arc_e = asg_get_idx_v_w_preedge_with_id(string_graph, x^1, ei, string_graph->ver[x^1].out_edge[wi[0]].adjvex, string_graph->ver[x^1].out_edge[wi[0]].id);
                    // assert(arc_e == 1);
                    // push_edge_list(&string_graph->ver[string_graph->ver[x^1].out_edge[wi[0]].adjvex].in_edge[ei[0]], vex, read_sta, read_end, edge_sta, edge_end);
                // }

                ret = reduce_edge_v_w(string_graph, v^1, w^1);
                assert(ret != -1);
                n_reduced++;
                break;
            }
        }
        if (identity < 0.9) // trim tip w->v
        {
            ret = reduce_edge_v_w(string_graph, w, v);
            assert(ret != -1);
            n_reduced++;
            ret = reduce_edge_v_w(string_graph, v^1, w^1);
            assert(ret != -1);
            n_reduced++;
        }
    }

    if (vi != NULL) free(vi);
    if (wi != NULL) free(wi);
    if (ei != NULL) free(ei);
    if (vseq != NULL)   free(vseq);
    if (xseq != NULL)   free(xseq);
    return n_reduced;
}

int olg_popping_small_bubble(StringGraph *string_graph)
{
    uint32_t v, *vi, w, *wi, *xi;
    uint32_t arc_n, arc_m_v_w, max_vern, arc_n_v_w;
    int edge_alloc = EDGE_TR * 2, n_bubble = 0, w_idx;
    uint32_t vertex_num = string_graph->nextId + 1;
    vi = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
    wi = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
    xi = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));

    for (v = 0; v < vertex_num; v++)
    {
        arc_n = get_v_arc_n(string_graph, v, vi);
        if (arc_n <= 1) continue;
        for (size_t i = 0; i < arc_n; i++)
        {
            max_vern = 0, w_idx = -1;
            w = string_graph->ver[v].out_edge[vi[i]].adjvex;
            arc_n_v_w = asg_get_idx_v_w_edge(string_graph, v, wi, w, 1);
            if (arc_n_v_w >= 2) // exsiting >2 v->w edge
            {
                // printf("%d:%d->%d\n", k++, v, w);
                n_bubble++;
                for (size_t j = 0; j < arc_n_v_w; j++)
                {
                    if (string_graph->ver[v].out_edge[wi[j]].edge_list.ver_n > max_vern)
                    {
                        max_vern = string_graph->ver[v].out_edge[wi[j]].edge_list.ver_n;
                        w_idx = j;
                    }
                }
                assert(w_idx != -1);
                for (size_t j = 0; j < arc_n_v_w; j++)
                {
                    if (j != w_idx)
                    {
                        string_graph->ver[v].out_edge[wi[j]].reduce = 1;
                    }
                }

                arc_m_v_w = asg_get_idx_v_w_preedge(string_graph, v, xi, w, 1);
                assert(arc_m_v_w == arc_n_v_w);
                for (size_t j = 0; j < arc_m_v_w; j++)
                {
                    if (string_graph->ver[w].in_edge[xi[j]].id != string_graph->ver[v].out_edge[wi[w_idx]].id)
                    {
                        string_graph->ver[w].in_edge[xi[j]].reduce = 1;
                    }
                }
            }
        }
    }

    if (vi != NULL) free(vi);
    if (wi != NULL) free(wi);
    if (xi != NULL) free(xi);

    return n_bubble;
}

int ver_in_cluster(uint32_t x, Cluster_t *cls)
{
    for (size_t i = 0; i < cls->n; i++)
    {
        if (x == cls->v[i])
            return 1;
    }
    return 0;
}

int ver_not_in_cluster(uint32_t x, Cluster_t *cls)
{
    for (size_t i = 0; i < cls->n; i++)
    {
        if (x == cls->v[i])
            return 0;
    }
    return 1;
}

int edge_not_in_path(uint32_t v, uint32_t w, Edge_List *path)
{
    for (size_t i = 0; i < path->n; i++)
    {
        if (v == path->v[i] && w == path->w[i])
            return 0;
    }
    return 1;
}

int finding_best_path(StringGraph *string_graph, uint32_t entry, uint32_t exit, Cluster_t *clus, uint32_t *path, Edge_List *edge_list)
{
    int edge_alloc = EDGE_TR * 2, k = -1;
    uint32_t *prev;
    uint32_t arc_n, path_n = 0;
    double *dist, min = 100.0, tmp_dist;
    uint32_t *vi, w, *mark;
    uint32_t vertex_num = (string_graph->nextId + 1)>>1;
    dist = (double *)calloc(vertex_num, sizeof(double));
    vi = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
    mark = (uint32_t *)calloc(vertex_num, sizeof(uint32_t));
    prev = (uint32_t *)calloc(vertex_num, sizeof(uint32_t));

    // init dist[]
    mark[entry>>1] = 1;
    for (size_t j = 0; j < vertex_num; j++)
    {
        dist[j] = 100.0;
        prev[j] = vertex_num + 1;
    }
    // entry -> w, update dist[w]
    arc_n = get_v_arc_n(string_graph, entry, vi);

    for (size_t j = 0; j < arc_n; j++)
    {
        w = string_graph->ver[entry].out_edge[vi[j]].adjvex;
        dist[w>>1] = 100.0/string_graph->ver[entry].out_edge[vi[j]].end_pos;
    }

    for (size_t i = 1; i < clus->n; i++)
    {
        min = 100.0;
        for (size_t j = 0; j < clus->n; j++)
        {
            if (mark[clus->v[j]>>1] == 0 && dist[clus->v[j]>>1] < min)
            {
                min = dist[clus->v[j]>>1];
                k = j;
            }
        }
        assert(k != -1);
        mark[clus->v[k]>>1] = 1;
        arc_n = get_v_arc_n(string_graph, clus->v[k], vi);

        for (size_t j = 0; j < arc_n; j++)
        {
            w = string_graph->ver[clus->v[k]].out_edge[vi[j]].adjvex;
            tmp_dist = min + 100.0 / string_graph->ver[clus->v[k]].out_edge[vi[j]].end_pos;
            if (mark[w>>1] == 0 && tmp_dist < dist[w>>1])
            {
                dist[w>>1] = tmp_dist;
                prev[w>>1] = clus->v[k];
            }
        }

    }

    // for (size_t j = 0; j < clus->n; j++)
    // {
    //     printf("dist[%d]=%f\n", clus->v[j], dist[clus->v[j]]);
    //     printf("prev[%d]=%d\n", clus->v[j], prev[clus->v[j]]);
    // }

    uint32_t cur_v = exit;
    while(cur_v != vertex_num + 1)
    {
        path[path_n++] = cur_v;
        cur_v = prev[cur_v>>1];
    }
    path[path_n++] = entry;

    for (int j = path_n-1; j > 0; j--)
    {
        edge_list->v[edge_list->n] = path[j];
        edge_list->w[edge_list->n++] = path[j-1];
        // printf("%d->%d\t", path[j].id, path[j-1].id);
    }
    // printf("\n");

    if (dist != NULL)   free(dist);
    if (vi != NULL)   free(vi);
    if (mark != NULL)   free(mark);
    if (prev != NULL)   free(prev);
    return path_n;
}

int olg_popping_super_bubble(StringGraph *string_graph, uint32_t len_tr)
{
    int n_bubble = 0, edge_alloc = EDGE_TR * 2, cluster_alloc = EDGE_TR * 2;
    uint32_t v, *vi, w, *wi, x, cur_v;
    uint32_t arc_n, arc_m, i, j, *mark, len;
    Cluster_t *cluster;
    Edge_List *edge, edge_list;
    int out_edge_inside, out_edge_outside, in_edge_inside, in_edge_outside;
    int ret, entry_n, exit_n, second_entry_n, second_exit_n;
    uint32_t entry = 0, exit = 0, *best_path, second_entry = 0, second_exit = 0;
    asgStack stk_v;
    uint32_t vertex_num = string_graph->nextId + 1;

    kv_init(stk_v);
    vi = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
    wi = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
    mark = (uint32_t *)calloc(vertex_num, sizeof(uint32_t));
    best_path = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
    cluster = (Cluster_t *)calloc(vertex_num, sizeof(Cluster_t));
    for (v = 0; v < vertex_num; v++)
    {
        cluster[v].n = 0;
        cluster[v].v = (uint32_t *)calloc(cluster_alloc, sizeof(uint32_t));
    }
    edge = (Edge_List *)calloc(vertex_num, sizeof(Edge_List));
    for (v = 0; v < vertex_num; v++)
    {
        edge[v].v = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
        edge[v].w = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
    }
    edge_list.v = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
    edge_list.w = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));

    // step 1: detect super bubble cluster
    for (v = 0; v < vertex_num; v++)
    {
        if (mark[v>>1] != 0)   continue;
        if (get_v_arc_n(string_graph, v, vi) == 0 && get_v_arc_m(string_graph, v, vi) == 0)   continue;
        edge_alloc = cluster_alloc = EDGE_TR * 2;
        kv_init(stk_v);
        kv_push(uint32_t, stk_v, v);
        mark[v>>1] = 1;
        cluster[v].v[cluster[v].n++] = v;
        // printf("%d\t", v);
        while (kv_size(stk_v) > 0)
        {
            cur_v = kv_pop(stk_v);
            arc_n = get_v_arc_n(string_graph, cur_v, vi);
            for (i = 0; i < arc_n; i++)
            {
                w = string_graph->ver[cur_v].out_edge[vi[i]].adjvex;
                len = string_graph->ver[cur_v].out_edge[vi[i]].end_pos;
                if (len < len_tr)
                {
                    edge[v].v[edge[v].n] = cur_v;
                    edge[v].w[edge[v].n++] = w;
                    if (edge[v].n == edge_alloc)
                    {
                        edge_alloc = edge_alloc << 1;
                        edge[v].v = (uint32_t *)realloc(edge[v].v, edge_alloc * sizeof(uint32_t));
                        edge[v].w = (uint32_t *)realloc(edge[v].w, edge_alloc * sizeof(uint32_t));
                    }
                    if (mark[w>>1] == 0)
                    {
                        kv_push(uint32_t, stk_v, w);
                        // printf("%d\t", w);
                        mark[w>>1] = 1;
                        cluster[v].v[cluster[v].n++] = w;
                        if (cluster[v].n == cluster_alloc)
                        {
                            cluster_alloc = cluster_alloc << 1;
                            cluster[v].v = (uint32_t*)realloc(cluster[v].v, cluster_alloc * sizeof(uint32_t));
                        }
                    }
                }
            }

            arc_m = get_v_arc_m(string_graph, cur_v, vi);
            for (i = 0; i < arc_m; i++)
            {
                w = string_graph->ver[cur_v].in_edge[vi[i]].adjvex;
                len = string_graph->ver[cur_v].in_edge[vi[i]].end_pos;
                if (len < len_tr)
                {
                    if (mark[w>>1] == 0)
                    {
                        kv_push(uint32_t, stk_v, w);
                        // printf("%d\t", w);
                        mark[w>>1] = 1;
                        cluster[v].v[cluster[v].n++] = w;
                        if (cluster[v].n == cluster_alloc)
                        {
                            cluster_alloc = cluster_alloc << 1;
                            cluster[v].v = (uint32_t*)realloc(cluster[v].v, cluster_alloc * sizeof(uint32_t));
                        }
                    }
                }
            }
        }
        // printf("\n");
    }

    // step 2: locate entry/exit vertices
    // entry vertex: in_edge_outside = 1, out_edge_inside = 1
    // exit vertex:  in_edge_inside = 1, out_edge_outside = 1
    for (v = 0; v < vertex_num; v++)
    {
        if (cluster[v].n > 3) // super bubble >= 4 vertices
        {
            // printf("%d\n",v>>1);
            // printf("cluster\t");
            // for (i = 0; i < cluster[v].n; i++)
            // {
            //     printf("%d\t",cluster[v].v[i]>>1);
            // }
            // printf("\n");

            // printf("edges\t");
            // for (i = 0; i < edge[v].n; i++)
            // {
            //     printf("%d->%d\t", edge[v].v[i]>>1, edge[v].w[i]>>1);
            // }
            // printf("\n");

            entry_n = exit_n = 0;
            second_entry_n = second_exit_n = 0;
            for (i = 0; i < cluster[v].n; i++)
            {
                in_edge_outside = out_edge_inside = 0;
                in_edge_inside = out_edge_outside = 0;
                w = cluster[v].v[i];
                // printf("%d\t", w>>1);

                arc_n = get_v_arc_n(string_graph, w, wi);
                // printf("%d out edge\n", arc_n);
                for (j = 0; j < arc_n; j++)
                {
                    x = string_graph->ver[w].out_edge[wi[j]].adjvex;
                    // printf("%d\n", x);
                    ret = ver_in_cluster(x, &cluster[v]);
                    if (ret == 1)   out_edge_inside = 1;
                    ret = ver_not_in_cluster(x, &cluster[v]);
                    if (ret == 1)   out_edge_outside = 1;
                    // printf("out_edge_inside=%d out_edge_outside=%d\n", out_edge_inside, out_edge_outside);
                }
                arc_m = get_v_arc_m(string_graph, w, wi);
                // printf("%d in edge\n", arc_m);
                for (j = 0; j < arc_m; j++)
                {
                    x = string_graph->ver[w].in_edge[wi[j]].adjvex;
                    // printf("%d %d\n", x, xt);
                    ret = ver_in_cluster(x, &cluster[v]);
                    if (ret == 1)   in_edge_inside = 1;
                    ret = ver_not_in_cluster(x, &cluster[v]);
                    if (ret == 1)   in_edge_outside = 1;
                    // printf("in_edge_inside=%d in_edge_outside=%d\n", in_edge_inside, in_edge_outside);
                }
 
                // printf("%d %d\t%d %d\n",in_edge_outside,out_edge_inside,in_edge_inside,out_edge_outside);
                if (in_edge_outside == 1 && out_edge_inside == 1) // entry
                {
                    entry_n++;
                    entry = w;
                }
                else if ((in_edge_outside == 0 && arc_m == 0) && out_edge_inside == 1) // entry of utg end
                {
                    second_entry_n++;
                    second_entry = w;
                }
                if (in_edge_inside == 1 && out_edge_outside == 1) // exit
                {
                    exit_n++;
                    exit = w;
                }
                else if ((out_edge_outside == 0 && arc_n == 0) && in_edge_inside == 1) // exit of utg end
                {
                    second_exit_n++;
                    second_exit = w;
                }
            }
            // printf("entry %d, exit %d\n", entry_n, exit_n);

            if (entry_n == 0 && second_entry_n == 1)
            {
                entry_n = second_entry_n;
                entry = second_entry;
            }
            if (exit_n == 0 && second_exit_n == 1)
            {
                exit_n = second_exit_n;
                exit = second_exit;
            }
            // printf("entry %d, exit %d\n", entry_n, exit_n);

            // step 3: find shortest path between entry and exit vertices
            if (entry_n == 1 && exit_n == 1)
            {
                edge_list.n = 0;
                // printf("entry %d, exit %d\n", entry>>1, exit>>1);
                ret = finding_best_path(string_graph, entry, exit, &cluster[v], best_path, &edge_list);
                if (ret != 0)   n_bubble++;
                // best_path stores id and type of the path, reduce other edge in the cluster
                // printf("retain edge\t");
                // for (i = 0; i < edge_list.n; i++)
                // {
                //     printf("%d->%d\t",edge_list.v[i]>>1,edge_list.w[i]>>1);
                // }
                // printf("\n");
                for (i = 0; i < edge[v].n; i++)
                {
                    ret = edge_not_in_path(edge[v].v[i], edge[v].w[i], &edge_list);
                    if (ret == 1) // not in path, reduce edge
                    {
                        // printf("reduce %d.%c->%d.%c\n", edge[v].v[i] >> 1, "BE"[edge[v].v[i] & 1], edge[v].w[i] >> 1, "BE"[edge[v].w[i] & 1]);
                        reduce_edge_v_w(string_graph, edge[v].v[i], edge[v].w[i]);
                        // printf("reduce %d.%c->%d.%c\n", (edge[v].w[i] ^ 1) >> 1, "BE"[(edge[v].w[i] ^ 1) & 1], (edge[v].v[i] ^ 1) >> 1, "BE"[(edge[v].v[i] ^ 1) & 1]);
                        reduce_edge_v_w(string_graph, edge[v].w[i]^1, edge[v].v[i]^1);
                    }
                }
            }
        }
    }

    kv_destroy(stk_v);
    if (vi != NULL) free(vi);
    if (wi != NULL) free(wi);
    if (mark != NULL) free(mark);
    if (best_path != NULL) free(best_path);

    for (v = 0; v < vertex_num; v++)
        if (cluster[v].v != NULL) free(cluster[v].v);
    if (cluster != NULL) free(cluster);

    if (edge_list.v != NULL) free(edge_list.v);
    if (edge_list.w != NULL) free(edge_list.w);
    for (v = 0; v < vertex_num; v++)
    {
        if (edge[v].v != NULL) free(edge[v].v);
        if (edge[v].w != NULL) free(edge[v].w);
    }
    if (edge != NULL) free(edge);
    return n_bubble;
}

// step 0: check unitig beg vertex, attach seq on it
uint32_t attach_unitig_beg_vertex(StringGraph *string_graph, READ_t *target)
{
    uint32_t edge_tr = EDGE_TR, n_base = 0;
    uint32_t v, w, id, *vi, *wi, *ei, arc_e, arc_n, seq_n;
    uint32_t vex, read_sta, read_end, edge_sta, edge_end, edge_len;
    uint32_t vertex_num = string_graph->nextId + 1;
    vi = (uint32_t *)calloc(edge_tr, sizeof(uint32_t));
    wi = (uint32_t *)calloc(edge_tr, sizeof(uint32_t));
    ei = (uint32_t *)calloc(edge_tr, sizeof(uint32_t));
    char *wseq, *xseq;
    wseq = (char *)calloc(MAX_SEQ_LEN_KSW, sizeof(char));
    xseq = (char *)calloc(MAX_SEQ_LEN_KSW, sizeof(char));

    for (v = 0; v < vertex_num; v++) // for each node v
    {
        arc_n = get_v_arc_n(string_graph, v, vi);
        // arc_m = get_v_arc_m(string_graph, v, wi);
        // if (arc_m != 0) continue; // v has 0 in edge
        if (arc_n >= 1) // if v has >=1 out edge, v->w
        {
            vex = v, read_sta = edge_sta = 0;
            read_end = edge_end = target[v>>1].read_length;
            edge_len = read_end - read_sta;
            for (size_t i = 0; i < arc_n; i++)
            {
                // printf("add seq %d %d->%d %d\n", v & 1, v >> 1, string_graph->ver[v].out_edge[vi[i]].adjvex & 1, string_graph->ver[v].out_edge[vi[i]].adjvex >> 1);
                string_graph->ver[v].out_edge[vi[i]].end_pos += edge_len;
                seq_n = cat_str(wseq, 0, string_graph->ver[v].out_edge[vi[i]].edge_seq, 0, string_graph->ver[v].out_edge[vi[i]].end_pos);
                if ((v&1)==0)
                    comp_cat_str(string_graph->ver[v].out_edge[vi[i]].edge_seq, 0, target[v>>1].read_seq, 0, read_end);
                else if ((v&1)==1)
                    cat_str(string_graph->ver[v].out_edge[vi[i]].edge_seq, 0, target[v>>1].read_seq, 0, read_end);
                cat_str(string_graph->ver[v].out_edge[vi[i]].edge_seq, read_end, wseq, 0, seq_n);
                n_base += edge_end;
                push_edge_list_at_beg(&string_graph->ver[v].out_edge[vi[i]], vex, read_sta, read_end, edge_sta, edge_end);
                arc_e = asg_get_idx_v_w_preedge_with_id(string_graph, v, ei, string_graph->ver[v].out_edge[vi[i]].adjvex, string_graph->ver[v].out_edge[vi[i]].id);
                assert(arc_e == 1);
                string_graph->ver[string_graph->ver[v].out_edge[vi[i]].adjvex].in_edge[ei[0]].end_pos += edge_len;
                push_edge_list_at_beg(&string_graph->ver[string_graph->ver[v].out_edge[vi[i]].adjvex].in_edge[ei[0]], vex, read_sta, read_end, edge_sta, edge_end);
            }
        }
    }

    for (v = 0; v < vertex_num; v++) // for each node v
    {
        arc_n = get_v_arc_n(string_graph, v, vi);
        if (arc_n >= 1) // if v has >=1 out edge, v->w
        {
            for (size_t i = 0; i < arc_n; i++)
            {
                if (string_graph->ver[v].out_edge[vi[i]].utg_id != 0)   continue;
                string_graph->ver[v].out_edge[vi[i]].utg_id = utg_id;
                string_graph->ver[v].out_edge[vi[i]].str = 0;
                w = string_graph->ver[v].out_edge[vi[i]].adjvex;
                id = string_graph->ver[v].out_edge[vi[i]].id;
                arc_e = asg_get_idx_v_w_preedge_with_id(string_graph, v, ei, w, id);
                printf("%d, %d %d -> %d %d:%d\n", arc_e, v & 1, v >> 1, w & 1, w >> 1, id);
                assert(arc_e == 1);
                string_graph->ver[w].in_edge[ei[0]].utg_id = utg_id;
                string_graph->ver[w].in_edge[ei[0]].str = 0;

                arc_e = asg_get_idx_v_w_edge(string_graph, w ^ 1, ei, v ^ 1, 1);
                printf("%d, %d %d -> %d %d\n", arc_e, (w ^ 1) & 1, (w ^ 1) >> 1, (v ^ 1) & 1, (v ^ 1) >> 1);
                assert(arc_e == 1);
                string_graph->ver[w^1].out_edge[ei[0]].utg_id = utg_id;
                string_graph->ver[w^1].out_edge[ei[0]].str = 1;
                id = string_graph->ver[w^1].out_edge[ei[0]].id;
                arc_e = asg_get_idx_v_w_preedge_with_id(string_graph, w ^ 1, ei, v ^ 1, id);
                assert(arc_e == 1);
                string_graph->ver[v^1].in_edge[ei[0]].utg_id = utg_id++;
                string_graph->ver[v^1].in_edge[ei[0]].str = 1;
            }
        }
    }

    if (vi != NULL) free(vi);
    if (wi != NULL) free(wi);
    if (ei != NULL) free(ei);
    if (wseq != NULL)   free(wseq);
    if (xseq != NULL)   free(xseq);
    return n_base;
}

int asg_cutting_branches(StringGraph *string_graph, StringGraph *bk_graph, READ_t *target)
{
    int ret, n_reduced = 0, end_D, end_I;
    double identity;
    uint32_t v, w, x, *vi, *wi, *ei, arc_n, arc_m, id, idx, newId;
    uint32_t edge_tr = EDGE_TR;
    uint32_t vertex_num = string_graph->nextId + 1;
    vi = (uint32_t *)calloc(edge_tr, sizeof(uint32_t));
    wi = (uint32_t *)calloc(edge_tr, sizeof(uint32_t));
    ei = (uint32_t *)calloc(edge_tr, sizeof(uint32_t));
    char *wseq, *xseq;
    wseq = (char *)calloc(MAX_SEQ_LEN_KSW, sizeof(char));
    xseq = (char *)calloc(MAX_SEQ_LEN_KSW, sizeof(char));
    uint32_t ksw_window = 2000;
    uint32_t vex, read_sta, read_end, edge_sta, edge_end;
    uint32_t seq_n, seq_sta_w, seq_sta_x, seqn_w, seqn_x, *seqn_cut, seqn_cut_value, cut_sta;
    strEdge *arc;
    arc = (strEdge *)calloc(edge_tr * 2, sizeof(strEdge));
    seqn_cut = (uint32_t *)calloc(edge_tr * 2, sizeof(uint32_t));
    int32_t in_idx, out_idx;

    // step 1: finding branches with out edge > 1
    for (v = 0; v < vertex_num; v++) // for each node v
    {
        arc_n = get_v_arc_n(string_graph, v, vi);
        arc_m = get_v_arc_m(string_graph, v, wi);
        if (arc_n > 1) // if v has > 1 out edge, v->w
        {
            for (size_t i = 0; i < arc_n; i++)
            {
                if (string_graph->ver[v].out_edge[vi[i]].edge_list.ver_n == 1)
                {
                    string_graph->ver[v].out_edge[vi[i]].reduce = 1;
                    w = string_graph->ver[v].out_edge[vi[i]].adjvex;
                    id = string_graph->ver[v].out_edge[vi[i]].id;
                    uint32_t arc_n2 = asg_get_idx_v_w_preedge_with_id(string_graph, v, ei, w, id);
                    assert(arc_n2 == 1);
                    string_graph->ver[w].in_edge[ei[0]].reduce = 1;

                    // TODO: if the ver_n of reverse edge != 1, copy its edge_list to other edges
                    arc_n2 = asg_get_idx_v_w_edge(string_graph, w ^ 1, ei, v ^ 1, 1);
                    assert(arc_n2 == 1);
                    string_graph->ver[w ^ 1].out_edge[ei[0]].reduce = 1;
                    id = string_graph->ver[w ^ 1].out_edge[ei[0]].id;
                    arc_n2 = asg_get_idx_v_w_preedge_with_id(string_graph, w ^ 1, ei, v ^ 1, id);
                    assert(arc_n2 == 1);
                    string_graph->ver[v ^ 1].in_edge[ei[0]].reduce = 1;
                }
            }
            arc_n = get_v_arc_n(string_graph, v, vi);
            if (arc_n <= 1) continue;

            printf("branching vertex %d\n", v >> 1);
            for (size_t i = 0; i < arc_n; i++)  arc[i] = string_graph->ver[v].out_edge[vi[i]];
            qsort(arc, arc_n, sizeof(strEdge), compare_arc);

            w = arc[0].adjvex;
            cut_sta = arc[0].edge_list.edge_sta[1];
            seqn_w = arc[0].end_pos - arc[0].edge_list.edge_end[0];
            printf("%d->%d:%d %d\n", v >> 1, w >> 1, seqn_w, arc[0].edge_list.ver_n);
            // finding cut position for each edge
            for (size_t i = 1; i < arc_n; i++) // edge v->w,v->x
            {
                x = arc[i].adjvex;
                seqn_x = arc[i].end_pos - arc[i].edge_list.edge_end[0];
                seqn_w = arc[0].end_pos - arc[0].edge_list.edge_end[0];
                printf("%d->%d:%d %d\n", v >> 1, x >> 1, seqn_x, arc[i].edge_list.ver_n);
                seq_n = ksw_window, seq_sta_w = arc[0].edge_list.edge_sta[1], seq_sta_x = arc[i].edge_list.edge_sta[1];
                seq_n = seq_n > seqn_x ? seqn_x : seq_n;
                cat_str(wseq, 0, arc[0].edge_seq, seq_sta_w, seq_n);
                cat_str(xseq, 0, arc[i].edge_seq, seq_sta_x, seq_n);
                identity = ksw2_extz_align(wseq, xseq, 1, -2, 2, 1, &end_D, &end_I);
                printf("%d identity: %f, %dD %dI\n", seq_n, identity, end_D, end_I);
                if (identity >= 0.9)
                {
                    while (seq_n == ksw_window)
                    {
                        seq_sta_w += seq_n, seq_sta_x += seq_n;
                        if (end_D > 0)  seq_sta_w -= end_D;
                        else if (end_I > 0) seq_sta_x -= end_I;
                        seqn_x = seqn_x - seq_n + end_I, seqn_w = seqn_w - seq_n + end_D;
                        seq_n = seq_n > seqn_x ? seqn_x : seq_n;
                        cat_str(wseq, 0, arc[0].edge_seq, seq_sta_w, seq_n);
                        cat_str(xseq, 0, arc[i].edge_seq, seq_sta_x, seq_n);
                        identity = ksw2_extz_align(wseq, xseq, 1, -2, 2, 1, &end_D, &end_I);
                        printf("%d identity: %f, %dD %dI\n", seq_n, identity, end_D, end_I);
                        if (identity < 0.9) break;
                    }
                }
                printf("%d %d %d %d\n",seq_sta_x, seq_sta_w, seqn_x, seqn_w);
                if (identity >= 0.9)    seqn_cut[i] = seq_sta_x + seqn_x;
                else    seqn_cut[i] = seq_sta_x;
            }

            // add the same sub sequence to v's in edge
            seqn_cut[0] = seqn_cut[1];
            qsort(seqn_cut, arc_n, sizeof(uint32_t), compare_uint32);
            seqn_cut_value = seqn_cut[0] - cut_sta;
            printf("cut_value = %d->%d\n", cut_sta + seqn_cut_value, seqn_cut_value);
            if (seqn_cut_value > 0)
            {
                if (arc_m == 0)
                {
                    newId = ++string_graph->nextId;
                    char *newSeq = (char *)calloc(MAX_SEQ_LEN_KSW, sizeof(char));
                    cat_str(newSeq, 0, arc[0].edge_seq, 0, cut_sta + seqn_cut_value);
                    asg_add_edge_v_w(string_graph, newId, v, 2, 0, cut_sta + seqn_cut_value, newSeq, &in_idx, &out_idx);
                    printf("1 add new edge %d->%d:%d\n", newId >> 1, v >> 1, cut_sta + seqn_cut_value);
                    string_graph->ver[newId].out_edge[out_idx].utg_id = utg_id;
                    string_graph->ver[newId].out_edge[out_idx].str = 0;
                    string_graph->ver[v].in_edge[in_idx].utg_id = utg_id;
                    string_graph->ver[v].in_edge[in_idx].str = 0;

                    for (size_t j = 0; j < arc[0].edge_list.ver_n; j++)
                    {
                        vex = arc[0].edge_list.vex[j];
                        read_sta = arc[0].edge_list.read_sta[j];
                        read_end = arc[0].edge_list.read_end[j];
                        edge_sta = arc[0].edge_list.edge_sta[j];
                        edge_end = arc[0].edge_list.edge_end[j];
                        if (edge_end <= seqn_cut_value + cut_sta)
                        {
                            push_edge_list(&string_graph->ver[newId].out_edge[out_idx], vex, read_sta, read_end, edge_sta, edge_end);
                            push_edge_list(&string_graph->ver[v].in_edge[in_idx], vex, read_sta, read_end, edge_sta, edge_end);
                        }
                        else if ((edge_sta < seqn_cut_value + cut_sta) && (edge_end > seqn_cut_value + cut_sta))
                        {
                            if ((vex&1)==0) read_sta = read_sta + (edge_end - (seqn_cut_value + cut_sta));
                            else if ((vex&1)==1)    read_end = read_end - (edge_end - (seqn_cut_value + cut_sta));
                            push_edge_list(&string_graph->ver[newId].out_edge[out_idx], vex, read_sta, read_end, edge_sta, seqn_cut_value + cut_sta);
                            push_edge_list(&string_graph->ver[v].in_edge[in_idx], vex, read_sta, read_end, edge_sta, seqn_cut_value + cut_sta);
                        }
                    }

                    newId = ++string_graph->nextId;
                    char *newSeqRev = (char *)calloc(MAX_SEQ_LEN_KSW, sizeof(char));
                    uint32_t arc_n2 = asg_get_idx_v_w_edge(string_graph, arc[0].adjvex ^ 1, vi, v ^ 1, 1);
                    assert(arc_n2 == 1);
                    int32_t cut_pos = string_graph->ver[arc[0].adjvex ^ 1].out_edge[vi[0]].end_pos - (cut_sta + seqn_cut_value);
                    uint32_t edge_len = cut_sta + seqn_cut_value;
                    if (cut_pos < 0)   cut_pos = 0, edge_len = string_graph->ver[arc[0].adjvex ^ 1].out_edge[vi[0]].end_pos; // the length of reverse edge is slightly different from forward edge
                    // printf("len = %d, cut_pos = %d, cut_len = %d\n", string_graph->ver[arc[0].adjvex ^ 1].out_edge[vi[0]].end_pos, cut_pos, edge_len);
                    cat_str(newSeqRev, 0, string_graph->ver[arc[0].adjvex ^ 1].out_edge[vi[0]].edge_seq, cut_pos, edge_len);
                    asg_add_edge_v_w(string_graph, v ^ 1, newId, 2, 0, cut_sta + seqn_cut_value, newSeqRev, &in_idx, &out_idx);
                    printf("1 add new edge %d->%d:%d\n", (v ^ 1) >> 1, newId >> 1, cut_sta + seqn_cut_value);
                    string_graph->ver[v ^ 1].out_edge[out_idx].utg_id = utg_id;
                    string_graph->ver[v ^ 1].out_edge[out_idx].str = 1;
                    string_graph->ver[newId].in_edge[in_idx].utg_id = utg_id++;
                    string_graph->ver[newId].in_edge[in_idx].str = 1;
                    for (size_t j = 0; j < string_graph->ver[arc[0].adjvex ^ 1].out_edge[vi[0]].edge_list.ver_n; j++)
                    {
                        vex = string_graph->ver[arc[0].adjvex ^ 1].out_edge[vi[0]].edge_list.vex[j];
                        read_sta = string_graph->ver[arc[0].adjvex ^ 1].out_edge[vi[0]].edge_list.read_sta[j];
                        read_end = string_graph->ver[arc[0].adjvex ^ 1].out_edge[vi[0]].edge_list.read_end[j];
                        edge_sta = string_graph->ver[arc[0].adjvex ^ 1].out_edge[vi[0]].edge_list.edge_sta[j];
                        edge_end = string_graph->ver[arc[0].adjvex ^ 1].out_edge[vi[0]].edge_list.edge_end[j];
                        if (edge_sta < cut_pos && edge_end >= cut_pos)
                        {
                            if ((vex&1)==1) read_sta = read_sta + (cut_pos - edge_sta);
                            else if ((vex&1)==0)    read_end = read_end - (cut_pos - edge_sta);
                            push_edge_list(&string_graph->ver[v ^ 1].out_edge[out_idx], vex, read_sta, read_end, 0, edge_end - cut_pos);
                            push_edge_list(&string_graph->ver[newId].in_edge[in_idx], vex, read_sta, read_end, 0, edge_end - cut_pos);
                        }
                        else if (edge_sta >=cut_pos)
                        {
                            push_edge_list(&string_graph->ver[v ^ 1].out_edge[out_idx], vex, read_sta, read_end, edge_sta - cut_pos, edge_end - cut_pos);
                            push_edge_list(&string_graph->ver[newId].in_edge[in_idx], vex, read_sta, read_end, edge_sta - cut_pos, edge_end - cut_pos);
                        }
                    }

                }
                else
                {
                    for (size_t i = 0; i < arc_m; i++) // edge w->v
                    {
                        uint32_t acc_len = string_graph->ver[v].in_edge[wi[i]].end_pos;
                        w = string_graph->ver[v].in_edge[wi[i]].adjvex;
                        id = string_graph->ver[v].in_edge[wi[i]].id;
                        cat_str(string_graph->ver[v].in_edge[wi[i]].edge_seq, string_graph->ver[v].in_edge[wi[i]].end_pos, arc[0].edge_seq, cut_sta, seqn_cut_value);
                        string_graph->ver[v].in_edge[wi[i]].end_pos += seqn_cut_value;
                        uint32_t arc_n2 = asg_get_idx_v_w_edge_with_id(string_graph, w, vi, v, id);
                        assert(arc_n2 == 1);
                        string_graph->ver[w].out_edge[vi[0]].end_pos += seqn_cut_value;
                        printf("2 add new seq to %d->%d:%d\n", w >> 1, v >> 1, seqn_cut_value);
                        for (size_t j = 0; j < arc[0].edge_list.ver_n; j++)
                        {
                            vex = arc[0].edge_list.vex[j];
                            read_sta = arc[0].edge_list.read_sta[j];
                            read_end = arc[0].edge_list.read_end[j];
                            edge_sta = arc[0].edge_list.edge_sta[j];
                            edge_end = arc[0].edge_list.edge_end[j];
                            if (edge_sta >= cut_sta && edge_end <= seqn_cut_value)
                            {
                                push_edge_list(&string_graph->ver[v].in_edge[wi[i]], vex, read_sta, read_end, edge_sta - cut_sta + acc_len, edge_end - cut_sta + acc_len);
                                push_edge_list(&string_graph->ver[w].out_edge[vi[0]], vex, read_sta, read_end, edge_sta - cut_sta + acc_len, edge_end - cut_sta + acc_len);
                            }
                            else if (edge_sta >= cut_sta && edge_sta < seqn_cut_value + cut_sta && edge_end > seqn_cut_value + cut_sta)
                            {
                                if ((vex & 1) == 0)
                                    read_sta = read_sta + (edge_end - (seqn_cut_value + cut_sta));
                                else if ((vex & 1) == 1)
                                    read_end = read_end - (edge_end - (seqn_cut_value + cut_sta));
                                push_edge_list(&string_graph->ver[v].in_edge[wi[i]], vex, read_sta, read_end, edge_sta - cut_sta + acc_len, seqn_cut_value + acc_len);
                                push_edge_list(&string_graph->ver[w].out_edge[vi[0]], vex, read_sta, read_end, edge_sta - cut_sta + acc_len, seqn_cut_value + acc_len);
                            }
                        }
                    }

                    uint32_t arc_n2 = asg_get_idx_v_w_edge(string_graph, arc[0].adjvex ^ 1, vi, v ^ 1, 1);
                    assert(arc_n2 == 1);
                    idx = vi[0];
                    uint32_t cut_pos = string_graph->ver[arc[0].adjvex ^ 1].out_edge[idx].end_pos - (cut_sta + seqn_cut_value);
                    for (size_t i = 0; i < arc_m; i++) // edge w->v
                    {
                        w = string_graph->ver[v].in_edge[wi[i]].adjvex;
                        arc_n2 = asg_get_idx_v_w_edge(string_graph, v ^ 1, vi, w ^ 1, 1);
                        assert(arc_n2 == 1);
                        cat_str(wseq, 0, string_graph->ver[v ^ 1].out_edge[vi[0]].edge_seq, 0, string_graph->ver[v ^ 1].out_edge[vi[0]].end_pos);
                        cat_str(string_graph->ver[v ^ 1].out_edge[vi[0]].edge_seq, 0, string_graph->ver[arc[0].adjvex ^ 1].out_edge[idx].edge_seq, cut_pos, seqn_cut_value);
                        cat_str(string_graph->ver[v ^ 1].out_edge[vi[0]].edge_seq, seqn_cut_value, wseq, 0, string_graph->ver[v ^ 1].out_edge[vi[0]].end_pos);
                        string_graph->ver[v ^ 1].out_edge[vi[0]].end_pos += seqn_cut_value;
                        arc_n2 = asg_get_idx_v_w_preedge_with_id(string_graph, v ^ 1, ei, w ^ 1, string_graph->ver[v ^ 1].out_edge[vi[0]].id);
                        assert(arc_n2 == 1);
                        string_graph->ver[w ^ 1].in_edge[ei[0]].end_pos += seqn_cut_value;
                        printf("2 add new seq to %d->%d:%d\n", (v ^ 1) >> 1, (w ^ 1) >> 1, seqn_cut_value);
                        for (int j = string_graph->ver[arc[0].adjvex ^ 1].out_edge[idx].edge_list.ver_n-1; j >= 0; j--)
                        {
                            vex = string_graph->ver[arc[0].adjvex ^ 1].out_edge[idx].edge_list.vex[j];
                            read_sta = string_graph->ver[arc[0].adjvex ^ 1].out_edge[idx].edge_list.read_sta[j];
                            read_end = string_graph->ver[arc[0].adjvex ^ 1].out_edge[idx].edge_list.read_end[j];
                            edge_sta = string_graph->ver[arc[0].adjvex ^ 1].out_edge[idx].edge_list.edge_sta[j];
                            edge_end = string_graph->ver[arc[0].adjvex ^ 1].out_edge[idx].edge_list.edge_end[j];
                            // printf("%d %d %d, %d %d\n", vex >> 1, edge_sta, edge_end, seqn_cut_value, cut_pos);
                            if (edge_sta >= cut_pos && edge_end <= cut_pos + seqn_cut_value)
                            {
                                push_edge_list_at_beg(&string_graph->ver[v ^ 1].out_edge[vi[0]], vex, read_sta, read_end, 0, edge_end - edge_sta);
                                push_edge_list_at_beg(&string_graph->ver[w ^ 1].in_edge[ei[0]], vex, read_sta, read_end, 0, edge_end - edge_sta);
                            }
                            else if (edge_sta >= cut_pos && edge_sta < cut_pos + seqn_cut_value && edge_end > cut_pos + seqn_cut_value)
                            {
                                if ((vex&1)==0) read_sta = read_sta + (edge_end - (cut_pos + seqn_cut_value));
                                else if ((vex&1)==1)    read_end = read_end - (edge_end - (cut_pos + seqn_cut_value));
                                push_edge_list_at_beg(&string_graph->ver[v ^ 1].out_edge[vi[0]], vex, read_sta, read_end, 0, cut_pos + seqn_cut_value - edge_sta);
                                push_edge_list_at_beg(&string_graph->ver[w ^ 1].in_edge[ei[0]], vex, read_sta, read_end, 0, cut_pos + seqn_cut_value - edge_sta);
                            }
                            else if (edge_sta < cut_pos && edge_end >= cut_pos && edge_end <= cut_pos + seqn_cut_value)
                            {
                                if ((vex&1)==1) read_sta = read_sta + (cut_pos - edge_sta);
                                else if ((vex&1)==0)    read_end = read_end - (cut_pos - edge_sta);
                                push_edge_list_at_beg(&string_graph->ver[v ^ 1].out_edge[vi[0]], vex, read_sta, read_end, 0, edge_end - cut_pos);
                                push_edge_list_at_beg(&string_graph->ver[w ^ 1].in_edge[ei[0]], vex, read_sta, read_end, 0, edge_end - cut_pos);
                            }
                            else if (edge_sta < cut_pos && edge_end >= cut_pos && edge_end > cut_pos + seqn_cut_value)
                            {
                                if ((vex&1)==1) read_sta = read_sta + (cut_pos - edge_sta), read_end = read_end - (edge_end - (cut_pos + seqn_cut_value));
                                else if ((vex&1)==0)    read_end = read_end - (cut_pos - edge_sta), read_sta = read_sta + (edge_end - (cut_pos + seqn_cut_value));
                                push_edge_list_at_beg(&string_graph->ver[v ^ 1].out_edge[vi[0]], vex, read_sta, read_end, 0, seqn_cut_value);
                                push_edge_list_at_beg(&string_graph->ver[w ^ 1].in_edge[ei[0]], vex, read_sta, read_end, 0, seqn_cut_value);
                            }
                        }
                    }
                }
            }

            // performing cutting process for each edge
            for (size_t i = 0; i < arc_n; i++) // edge v->w,v->x
            {
                x = arc[i].adjvex;
                ret = reduce_edge_v_w(string_graph, v, x);
                if (ret != -1)  n_reduced++;
                printf("1 reduce edge %d->%d\n", v >> 1, x >> 1);

                if (abs(seqn_cut_value + cut_sta - arc[i].end_pos) <= 100)
                {
                    ret = reduce_edge_v_w(string_graph, x ^ 1, v ^ 1);
                    if (ret != -1)  n_reduced++;
                    printf("2 reduce edge %d->%d\n", x >> 1, v >> 1);
                    continue;
                }

                newId = ++string_graph->nextId;
                char *newSeq = (char *)calloc(MAX_SEQ_LEN_KSW, sizeof(char));
                cat_str(newSeq, 0, arc[i].edge_seq, cut_sta + seqn_cut_value, arc[i].end_pos);
                asg_add_edge_v_w(string_graph, newId, x, 2, 0, arc[i].end_pos - (seqn_cut_value + cut_sta), newSeq, &in_idx, &out_idx);
                printf("3 add new edge %d->%d: %d\n",newId>>1,x>>1,arc[i].end_pos - (seqn_cut_value + cut_sta));
                string_graph->ver[newId].out_edge[out_idx].utg_id = arc[i].utg_id;
                string_graph->ver[newId].out_edge[out_idx].str = arc[i].str;
                string_graph->ver[x].in_edge[in_idx].utg_id = arc[i].utg_id;
                string_graph->ver[x].in_edge[in_idx].str = arc[i].str;
                for (size_t j = 0; j < arc[i].edge_list.ver_n; j++)
                {
                    vex = arc[i].edge_list.vex[j];
                    read_sta = arc[i].edge_list.read_sta[j];
                    read_end = arc[i].edge_list.read_end[j];
                    edge_sta = arc[i].edge_list.edge_sta[j];
                    edge_end = arc[i].edge_list.edge_end[j];
                    if (edge_sta < cut_sta + seqn_cut_value && edge_end > cut_sta + seqn_cut_value)
                    {
                        if ((vex&1)==1) read_sta = read_sta + (cut_sta + seqn_cut_value - edge_sta);
                        else if ((vex&1)==0)    read_end = read_end - (cut_sta + seqn_cut_value - edge_sta);
                        push_edge_list(&string_graph->ver[newId].out_edge[out_idx], vex, read_sta, read_end, 0, edge_end - (cut_sta + seqn_cut_value));
                        push_edge_list(&string_graph->ver[x].in_edge[in_idx], vex, read_sta, read_end, 0, edge_end - (cut_sta + seqn_cut_value));
                    }
                    else if (edge_sta >= cut_sta + seqn_cut_value)
                    {
                        edge_sta -= (cut_sta + seqn_cut_value);
                        edge_end -= (cut_sta + seqn_cut_value);
                        push_edge_list(&string_graph->ver[newId].out_edge[out_idx], vex, read_sta, read_end, edge_sta, edge_end);
                        push_edge_list(&string_graph->ver[x].in_edge[in_idx], vex, read_sta, read_end, edge_sta, edge_end);
                    }
                }

                uint32_t arc_n2 = asg_get_idx_v_w_edge(string_graph, x ^ 1, vi, v ^ 1, 1);
                assert(arc_n2 == 1);
                idx = vi[0];
                uint32_t cut_pos = string_graph->ver[x ^ 1].out_edge[idx].end_pos - (cut_sta + seqn_cut_value);

                newId = ++string_graph->nextId;
                char *newSeqRev = (char *)calloc(MAX_SEQ_LEN_KSW, sizeof(char));
                cat_str(newSeqRev, 0, string_graph->ver[x ^ 1].out_edge[idx].edge_seq, 0, cut_pos);
                asg_add_edge_v_w(string_graph, x ^ 1, newId, 2, 0, cut_pos, newSeqRev, &in_idx, &out_idx);
                printf("3 add new edge %d->%d: %d\n", x >> 1, newId >> 1, cut_pos);
                string_graph->ver[x ^ 1].out_edge[out_idx].utg_id = string_graph->ver[x ^ 1].out_edge[idx].utg_id;
                string_graph->ver[x ^ 1].out_edge[out_idx].str = string_graph->ver[x ^ 1].out_edge[idx].str;
                string_graph->ver[newId].in_edge[in_idx].utg_id = string_graph->ver[x ^ 1].out_edge[idx].utg_id;
                string_graph->ver[newId].in_edge[in_idx].str = string_graph->ver[x ^ 1].out_edge[idx].str;
                for (size_t j = 0; j < string_graph->ver[x ^ 1].out_edge[idx].edge_list.ver_n; j++)
                {
                    vex = string_graph->ver[x ^ 1].out_edge[idx].edge_list.vex[j];
                    read_sta = string_graph->ver[x ^ 1].out_edge[idx].edge_list.read_sta[j];
                    read_end = string_graph->ver[x ^ 1].out_edge[idx].edge_list.read_end[j];
                    edge_sta = string_graph->ver[x ^ 1].out_edge[idx].edge_list.edge_sta[j];
                    edge_end = string_graph->ver[x ^ 1].out_edge[idx].edge_list.edge_end[j];
                    if (edge_end <= cut_pos)
                    {
                        push_edge_list(&string_graph->ver[x ^ 1].out_edge[out_idx], vex, read_sta, read_end, edge_sta, edge_end);
                        push_edge_list(&string_graph->ver[newId].in_edge[in_idx], vex, read_sta, read_end, edge_sta, edge_end);
                    }
                    else if (edge_sta < cut_pos && edge_end > cut_pos)
                    {
                        if ((vex&1)==0) read_sta = read_sta + (edge_end - cut_pos);
                        else if ((vex&1)==1) read_end = read_end - (edge_end - cut_pos);
                        push_edge_list(&string_graph->ver[x ^ 1].out_edge[out_idx], vex, read_sta, read_end, edge_sta, cut_pos);
                        push_edge_list(&string_graph->ver[newId].in_edge[in_idx], vex, read_sta, read_end, edge_sta, cut_pos);
                    }
                }
                ret = reduce_edge_v_w(string_graph, x ^ 1, v ^ 1);
                if (ret != -1)  n_reduced++;
                printf("3 reduce edge %d->%d\n", x >> 1, v >> 1);
            }
        }
    }

    if (arc != NULL)  free(arc);
    if (seqn_cut != NULL)  free(seqn_cut);

    if (vi != NULL) free(vi);
    if (wi != NULL) free(wi);
    if (ei != NULL) free(ei);
    if (wseq != NULL)   free(wseq);
    if (xseq != NULL)   free(xseq);
    return n_reduced;
}

uint32_t output_unitig(UTG_t *utg, const char *out_file_dir, uint32_t n_utg)
{
    char name[32];
    FILE *fp_res;

    fp_res = fopen(out_file_dir, "w");
    if (fp_res == NULL)
    {
        fprintf(stderr, "[Wrong] Failed to open file %s!!!\n", out_file_dir);
        exit(1);
    }

    for (int32_t i = 0; i < n_utg; i++)
    {
        if (utg[i].len > 0)
        {
            sprintf(name, ">utg%c%d", "+-"[utg[i].str], utg[i].id);
            fprintf(fp_res, "%s\n", name);
            fprintf(fp_res, "%s\n", utg[i].edge_seq);
        }
    }

    fclose(fp_res);
    return 0;
}

uint32_t in_utg(uint32_t vex, uint32_t vi, uint32_t ei,StringGraph *string_graph)
{
    uint32_t bflag = 0, arc_n;
    uint32_t vertex_num = string_graph->nextId + 1;
    for (int v = 0; v < vertex_num; v++)
    {
        arc_n = string_graph->ver[v].arc_n;
        for (int i = 0; i < arc_n; i++)
        {
            if (string_graph->ver[v].out_edge[i].reduce != 1 && string_graph->ver[v].out_edge[i].edge_list.ver_n > 1)
            {
                for (int k = 0; k < string_graph->ver[v].out_edge[i].edge_list.ver_n; k++)
                {
                    if (string_graph->ver[v].out_edge[i].edge_list.vex[k]>>1 == vex>>1 && (!(v == vi && i == ei)))
                    {
                        string_graph->ver[vi].out_edge[ei].reduce = 1;
                        bflag = 1;
                        break;
                    }
                }
            }
            if (bflag == 1) break;
        }
        if (bflag == 1) break;
    }
    return bflag;
}

int generating_unitig(StringGraph *string_graph, UTG_t *utg, uint32_t tr, uint32_t vern_tr)
{
    uint32_t w, id, idx, n_utg = 0, arc_n = 0;
    uint32_t vertex_num = string_graph->nextId + 1;
    uint32_t *vi = (uint32_t *)calloc(EDGE_TR, sizeof(uint32_t));

    for (uint32_t v = 0; v < vertex_num; v++)
    {
        for (uint32_t i = 0; i < string_graph->ver[v].arc_n; i++)
        {
            if (string_graph->ver[v].out_edge[i].reduce == 1)   continue;
            // filtering utg containing 1 read which contained in other longer utg
            // if ((string_graph->ver[v].out_edge[i].edge_list.ver_n <= vern_tr) || string_graph->ver[v].out_edge[i].end_pos < tr)
            if ((string_graph->ver[v].out_edge[i].edge_list.ver_n <= vern_tr && in_utg(string_graph->ver[v].out_edge[i].edge_list.vex[0], v, i, string_graph) == 1) || string_graph->ver[v].out_edge[i].end_pos < tr)
            {
                w = string_graph->ver[v].out_edge[i].adjvex;
                id = string_graph->ver[v].out_edge[i].id;
                string_graph->ver[v].out_edge[i].reduce = 1;
                arc_n = asg_get_idx_v_w_preedge_with_id(string_graph, v, vi, w, id);
                assert(arc_n == 1);
                string_graph->ver[w].in_edge[vi[0]].reduce = 1;
                arc_n = asg_get_idx_v_w_edge(string_graph, w ^ 1, vi, v ^ 1, 1);
                assert(arc_n == 1);
                string_graph->ver[w ^ 1].out_edge[vi[0]].reduce = 1;
                id = string_graph->ver[w ^ 1].out_edge[vi[0]].id;
                arc_n = asg_get_idx_v_w_preedge_with_id(string_graph, w ^ 1, vi, v ^ 1, id);
                assert(arc_n == 1);
                string_graph->ver[v ^ 1].in_edge[vi[0]].reduce = 1;
                printf("reduce edge %d.%c->%d.%c\n", v >> 1, "BE"[v & 1], w >> 1, "BE"[w & 1]);
                printf("reduce edge %d.%c->%d.%c\n", (w ^ 1) >> 1, "BE"[(w ^ 1) & 1], (v ^ 1) >> 1, "BE"[(v ^ 1) & 1]);
            }
        }
    }

    for (uint32_t v = 0; v < vertex_num; v++)
    {
        for (uint32_t i = 0; i < string_graph->ver[v].arc_n; i++)
        {
            if (string_graph->ver[v].out_edge[i].reduce != 1)
            {
                idx = (string_graph->ver[v].out_edge[i].utg_id << 1) ^ string_graph->ver[v].out_edge[i].str;
                utg[idx].id = string_graph->ver[v].out_edge[i].utg_id;
                utg[idx].str = string_graph->ver[v].out_edge[i].str;
                utg[idx].len = string_graph->ver[v].out_edge[i].end_pos;
                utg[idx].seq_n = string_graph->ver[v].out_edge[i].end_pos;
                utg[idx].edge_list = string_graph->ver[v].out_edge[i].edge_list;
                if (string_graph->ver[v].out_edge[i].utg_id > n_utg)    n_utg = string_graph->ver[v].out_edge[i].utg_id;
                // utg[idx].edge_seq = string_graph->ver[v].out_edge[i].edge_seq;
                char *seq;
                seq = (char *)calloc(utg[idx].seq_n + 1, sizeof(char));
                cat_str(seq, 0, string_graph->ver[v].out_edge[i].edge_seq, 0, string_graph->ver[v].out_edge[i].end_pos);
                utg[idx].edge_seq = seq;
                printf("utg %c%d, len %d %ld, ver_n %d\n", "+-"[utg[idx].str], utg[idx].id, utg[idx].len, strlen(utg[idx].edge_seq), utg[idx].edge_list.ver_n);
            }
        }
    }

    if (vi != NULL) free(vi);
    return (n_utg + 1) << 1;
}

int create_assembly_graph(StringGraph *ol_string_graph, StringGraph *bk_string_graph, OVE_C_t *ove_cl, READ_t *target, UTG_t *utg, const char *out_file_dir)
{
    int iter = 5;
    uint32_t i, j, f, g, ft, gt, tmp, sta_pos, end_pos, step, w, arc_n;
    uint32_t tl3 = 0, ql3 = 0;
    uint32_t tl5 = 0, ql5 = 0;
    uint32_t n_skip = 0, n_edge = 0, n_reduced = 0, n_merged = 0, n_bubble = 0, n_utg = 0;
    uint32_t edge_tr = EDGE_TR, edge_alloc = edge_tr * 2; // TODO: finding appropriate value
    uint32_t tip_tr = 10000, super_bubble_tr = 10000; // if tip > 10000 bp, not trim to retain information
    uint32_t *vi = (uint32_t *)calloc(edge_alloc, sizeof(uint32_t));
    uint32_t vertex_num = ol_string_graph->allocN;
    int32_t in_idx, out_idx;

    for (i = 0; i < ove_cl->n; ++i)
    {
        ql3 = ove_cl->ove[i].qs;
        tl3 = ove_cl->ove[i].rev == 0 ? ove_cl->ove[i].ts : ove_cl->ove[i].tl - ove_cl->ove[i].ts + 1;
        ql5 = ove_cl->ove[i].ql - ove_cl->ove[i].qe;
        tl5 = ove_cl->ove[i].rev == 0 ? ove_cl->ove[i].tl - ove_cl->ove[i].te : ove_cl->ove[i].te;
        // if (ove_cl->ove[i].qid == 803 || ove_cl->ove[i].tid == 803)
        // {
        //     printf("%d->%d\n", ove_cl->ove[i].tid, ove_cl->ove[i].qid);
        //     printf("%d\t%d\t%d\t%d\n", ove_cl->ove[i].ts, ove_cl->ove[i].te, ove_cl->ove[i].qs, ove_cl->ove[i].qe);
        // }
        // skip contained reads
        if ((ql3 <= tl3 && ql5 <= tl5) || (tl3 <= ql3 && tl5 <= ql5)) {n_skip++; continue;}
        f = ove_cl->ove[i].qid;
        g = ove_cl->ove[i].tid;
        if (ove_cl->ove[i].te < ove_cl->ove[i].ts)
        {
            tmp = ove_cl->ove[i].te;
            ove_cl->ove[i].te = ove_cl->ove[i].ts;
            ove_cl->ove[i].ts = tmp;
        }

        if (ql3 > tl3)
        {
            //      f.B        f.E
            //   f  ----------->
            //   g         ------------->
            //             g.B          g.E
            if (ove_cl->ove[i].rev == 0)
            {
                // add edge f.E -> g.E, f.E-> out edge, ->g.E in edge
                n_edge++;
                ft = 1, gt = 1;
                sta_pos = ove_cl->ove[i].te + ove_cl->ove[i].ql - ove_cl->ove[i].qe;
                end_pos = ove_cl->ove[i].tl;
                if (ol_string_graph->ver[(f << 1) ^ ft].arc_n >= edge_tr || ol_string_graph->ver[(g << 1) ^ gt].arc_m >= edge_tr)   continue;
                asg_add_edge_v_w(ol_string_graph, (f << 1) ^ ft, (g << 1) ^ gt, 0, sta_pos, end_pos, NULL, &in_idx, &out_idx);
                // printf("%d.E->%d.E:[%d-%d] %d\n", f, g, sta_pos, end_pos, end_pos - sta_pos);

                // add edge g.B -> f.B, g.B-> out edge, ->f.B in edge
                n_edge++;
                ft = 0, gt = 0;
                sta_pos = 0;
                end_pos = ove_cl->ove[i].qs - ove_cl->ove[i].ts;
                asg_add_edge_v_w(ol_string_graph, (g << 1) ^ gt, (f << 1) ^ ft, 0, sta_pos, end_pos, NULL, &in_idx, &out_idx);
                // printf("%d.B->%d.B:[%d-%d] %d\n", g, f, sta_pos, end_pos, end_pos - sta_pos);
            }
            //      f.B        f.E
            //   f  ----------->
            //   g         <-------------
            //             g.E          g.B
            else if (ove_cl->ove[i].rev == 1)
            {
                // add edge f.E -> g.B, f.E-> out edge, ->g.B in edge
                n_edge++;
                ft = 1, gt = 0;
                sta_pos = 0;
                end_pos = ove_cl->ove[i].ts - (ove_cl->ove[i].ql - ove_cl->ove[i].qe);
                if (ol_string_graph->ver[(f << 1) ^ ft].arc_n >= edge_tr || ol_string_graph->ver[(g << 1) ^ gt].arc_m >= edge_tr)   continue;
                asg_add_edge_v_w(ol_string_graph, (f << 1) ^ ft, (g << 1) ^ gt, 0, sta_pos, end_pos, NULL, &in_idx, &out_idx);
                // printf("%d.E->%d.B:[%d-%d] %d\n", f, g, sta_pos, end_pos, end_pos - sta_pos);

                // add edge g.E -> f.B, g.E-> out edge, ->f.B in edge
                n_edge++;
                ft = 0, gt = 1;
                sta_pos = 0;
                end_pos = ove_cl->ove[i].qs - (ove_cl->ove[i].tl - ove_cl->ove[i].te);
                asg_add_edge_v_w(ol_string_graph, (g << 1) ^ gt, (f << 1) ^ ft, 0, sta_pos, end_pos, NULL, &in_idx, &out_idx);
                // printf("%d.E->%d.B:[%d-%d] %d\n", g, f, sta_pos, end_pos, end_pos - sta_pos);
            }
        }
        else if (ql3 < tl3)
        {
            //               f.B        f.E
            //   f           ----------->
            //   g   ------------->
            //       g.B          g.E
            if (ove_cl->ove[i].rev == 0)
            {
                // add edge f.B -> g.B, f.B-> out edge, ->g.B in edge
                n_edge++;
                ft = 0, gt = 0;
                sta_pos = 0;
                end_pos = ove_cl->ove[i].ts - ove_cl->ove[i].qs;
                if (ol_string_graph->ver[(f << 1) ^ ft].arc_n >= edge_tr || ol_string_graph->ver[(g << 1) ^ gt].arc_m >= edge_tr)   continue;
                asg_add_edge_v_w(ol_string_graph, (f << 1) ^ ft, (g << 1) ^ gt, 0, sta_pos, end_pos, NULL, &in_idx, &out_idx);
                // printf("%d.B->%d.B:[%d-%d] %d\n", f, g, sta_pos, end_pos, end_pos - sta_pos);

                // add edge g.E -> f.E, g.E-> out edge, ->f.E in edge
                n_edge++;
                ft = 1, gt = 1;
                sta_pos = ove_cl->ove[i].qe + ove_cl->ove[i].tl - ove_cl->ove[i].te;
                end_pos = ove_cl->ove[i].ql;
                asg_add_edge_v_w(ol_string_graph, (g << 1) ^ gt, (f << 1) ^ ft, 0, sta_pos, end_pos, NULL, &in_idx, &out_idx);
                // printf("%d.E->%d.E:[%d-%d] %d\n", g, f, sta_pos, end_pos, end_pos - sta_pos);
            }
            //                 f.B        f.E
            //   f             ----------->
            //   g     <-------------
            //         g.E          g.B
            else if (ove_cl->ove[i].rev == 1)
            {
                // add edge f.B -> g.E, f.B-> out edge, ->g.E in edge
                n_edge++;
                ft = 0, gt = 1;
                sta_pos = ove_cl->ove[i].te + ove_cl->ove[i].qs;
                end_pos = ove_cl->ove[i].tl;
                if (ol_string_graph->ver[(f << 1) ^ ft].arc_n >= edge_tr || ol_string_graph->ver[(g << 1) ^ gt].arc_m >= edge_tr)   continue;
                asg_add_edge_v_w(ol_string_graph, (f << 1) ^ ft, (g << 1) ^ gt, 0, sta_pos, end_pos, NULL, &in_idx, &out_idx);
                // printf("%d.B->%d.E:[%d-%d] %d\n", f, g, sta_pos, end_pos, end_pos - sta_pos);

                // add edge g.B -> f.E, g.B-> out edge, ->f.E in edge
                n_edge++;
                ft = 1, gt = 0;
                sta_pos = ove_cl->ove[i].qe + ove_cl->ove[i].ts;
                end_pos = ove_cl->ove[i].ql;
                asg_add_edge_v_w(ol_string_graph, (g << 1) ^ gt, (f << 1) ^ ft, 0, sta_pos, end_pos, NULL, &in_idx, &out_idx);
                // printf("%d.B->%d.E:[%d-%d] %d\n", g, f, sta_pos, end_pos, end_pos - sta_pos);
            }
        }
    }
    uint32_t edge_sta = 0;
    for (i = 0; i < vertex_num; ++i)
    {
        for (j = 0; j < ol_string_graph->ver[i].arc_n; j++)
        {
            edge_sta = 0;
            push_edge_list(&ol_string_graph->ver[i].out_edge[j], ol_string_graph->ver[i].out_edge[j].adjvex, ol_string_graph->ver[i].out_edge[j].sta_pos, ol_string_graph->ver[i].out_edge[j].end_pos, edge_sta, edge_sta+ol_string_graph->ver[i].out_edge[j].end_pos - ol_string_graph->ver[i].out_edge[j].sta_pos);
            w = ol_string_graph->ver[i].out_edge[j].adjvex;
            arc_n = asg_get_idx_v_w_preedge_with_id(ol_string_graph, i, vi, w, ol_string_graph->ver[i].out_edge[j].id);
            assert(arc_n == 1);
            push_edge_list(&ol_string_graph->ver[w].in_edge[vi[0]], ol_string_graph->ver[i].out_edge[j].adjvex, ol_string_graph->ver[i].out_edge[j].sta_pos, ol_string_graph->ver[i].out_edge[j].end_pos, edge_sta, edge_sta+ol_string_graph->ver[i].out_edge[j].end_pos - ol_string_graph->ver[i].out_edge[j].sta_pos);
        }
        if (ol_string_graph->ver[i].arc_n > edge_tr)    printf("%d n %d\n",i>>1,ol_string_graph->ver[i].arc_n);
    }
    fprintf(stderr, "[Graph Construction] total %d edges\n", n_edge);
    fprintf(stderr, "[Graph Construction] skipping %d contained reads\n", n_skip);

    // step 0: generating string graph with breakpoint vertex: transitive -> merge path
    step = 0;
    n_reduced = olg_transitive_reduction(ol_string_graph, step);
    if (n_reduced != 0) fprintf(stderr, "[Graph Cleaning] transitively reduced %d arcs\n", n_reduced);

    n_merged = olg_merge_unbranching_path(ol_string_graph, target, step);
    if (n_merged != 0) fprintf(stderr, "[Path Merging] merged %d arcs in un-branching path\n", n_merged);

    // step 1: cleaning string graph for first iteration: tip -> transitive
    step = 1, n_reduced = 1, iter = 5;
    while (n_reduced != 0 && iter-- > 0)
    {
        n_reduced = olg_trimming_tips(ol_string_graph, tip_tr);
        if (n_reduced != 0) fprintf(stderr, "[Graph Cleaning] merged %d tips\n", n_reduced);

        n_merged = olg_merge_unbranching_path(ol_string_graph, target, step);
        if (n_merged != 0) fprintf(stderr, "[Path Merging] merged %d arcs in un-branching path\n", n_merged);

        n_merged += olg_transitive_reduction(ol_string_graph, step);
        if (n_merged != 0) fprintf(stderr, "[Graph Cleaning] transitively reduced %d arcs\n", n_merged);
    }

    // step 2: cleaning string graph for first iteration: merge path -> small bubble -> nested bubble
    n_bubble = 1;
    while (n_bubble != 0)
    {
        n_merged = olg_merge_unbranching_path(ol_string_graph, target, step);
        if (n_merged != 0) fprintf(stderr, "[Path Merging] merged %d arcs in un-branching path\n", n_merged);

        n_bubble = olg_popping_small_bubble(ol_string_graph);
        if (n_bubble != 0) fprintf(stderr, "[Graph Cleaning] merged %d bubbles\n", n_bubble);
    }

    // step 3: cleaning string graph for second iteration: merge path -> super bubble
    i = 1, n_bubble = 1, iter = 5;
    while (n_bubble != 0 && iter-- > 0)
    {
        n_merged = olg_merge_unbranching_path(ol_string_graph, target, step);
        if (n_merged != 0) fprintf(stderr, "[Path Merging] merged %d arcs in un-branching path\n", n_merged);

        n_bubble = olg_popping_super_bubble(ol_string_graph, super_bubble_tr*i); i++;
        if (n_bubble != 0) fprintf(stderr, "[Graph Cleaning] merged %d super bubbles\n", n_bubble);

        n_reduced = olg_trimming_tips(ol_string_graph, tip_tr);
        if (n_reduced != 0) fprintf(stderr, "[Graph Cleaning] merged %d tips\n", n_reduced);

        n_reduced = olg_popping_small_bubble(ol_string_graph);
        if (n_reduced != 0) fprintf(stderr, "[Graph Cleaning] merged %d bubbles\n", n_reduced);
    }

    // step 4: deal with the potentially wrong connection: out branches -> in branches
    n_merged = olg_merge_unbranching_path(ol_string_graph, target, step);
    if (n_merged != 0) fprintf(stderr, "[Path Merging] merged %d arcs in un-branching path\n", n_merged);

    uint32_t n_base = attach_unitig_beg_vertex(ol_string_graph, target);
    fprintf(stderr, "[Graph Cleaning] attaching %d sequence in unitig beg vertex\n", n_base);

    n_reduced = asg_cutting_branches(ol_string_graph, bk_string_graph, target);
    if (n_reduced != 0) fprintf(stderr, "[Graph Cleaning] cutting %d connect arcs\n", n_reduced);

    // step 5: finally merge path
    step = 2;
    n_merged = olg_merge_unbranching_path(ol_string_graph, target, step);
    if (n_merged != 0) fprintf(stderr, "[Path Merging] merged %d arcs in un-branching path\n", n_merged);

    n_utg = generating_unitig(ol_string_graph, utg, 1000, 1);

    output_unitig(utg, out_file_dir, n_utg);
    if (n_utg != 0) fprintf(stderr, "[Output] output %d unitigs in file\n", n_utg);

#ifdef _DRAW
    FILE *fp_g;
    char temp_g_dir[1024];
    memset(temp_g_dir, 0, 1024);
    strcpy(temp_g_dir, "./graph.txt");
    fp_g = fopen(temp_g_dir, "w");
    if (fp_g == NULL)
    {
        fprintf(stderr, "[%s] Failed to open file %s!!!\n", __func__, temp_g_dir);
        exit(1);
    }

    for (i = 0; i <= ol_string_graph->nextId; i++)
    {
        for (int j = 0; j < ol_string_graph->ver[i].arc_n; j++)
        {
            if (ol_string_graph->ver[i].out_edge[j].reduce != 1)
            {
                fprintf(fp_g, "%d.%c,", i >> 1, "BE"[i & 1]);
                fprintf(fp_g, "%d.%c,", ol_string_graph->ver[i].out_edge[j].adjvex >> 1, "BE"[ol_string_graph->ver[i].out_edge[j].adjvex & 1]);
                fprintf(fp_g, "%d,%d,%d,", ol_string_graph->ver[i].out_edge[j].end_pos, ol_string_graph->ver[i].out_edge[j].reduce, ol_string_graph->ver[i].out_edge[j].id);
                fprintf(fp_g, "%c%d,%d,", "+-"[ol_string_graph->ver[i].out_edge[j].str], ol_string_graph -> ver[i].out_edge[j].utg_id, ol_string_graph->ver[i].out_edge[j].edge_list.ver_n);
                for (size_t k = 0; k < ol_string_graph->ver[i].out_edge[j].edge_list.ver_n; k++)
                {
                    fprintf(fp_g, "[%d.%c %d-%d,%d-%d] ", ol_string_graph->ver[i].out_edge[j].edge_list.vex[k] >> 1, "BE"[ol_string_graph->ver[i].out_edge[j].edge_list.vex[k] & 1], ol_string_graph -> ver[i].out_edge[j].edge_list.edge_sta[k], ol_string_graph->ver[i].out_edge[j].edge_list.edge_end[k], ol_string_graph->ver[i].out_edge[j].edge_list.read_sta[k], ol_string_graph->ver[i].out_edge[j].edge_list.read_end[k]);
                }
                fprintf(fp_g, "\n");
            }
        }
    }

    // fprintf(fp_g, "\n");
    // for (i = 0; i <= ol_string_graph->nextId; i++)
    // {
    //     for (int j = 0; j < ol_string_graph->ver[i].arc_m; j++)
    //     {
    //         if (ol_string_graph->ver[i].in_edge[j].reduce != 1)
    //         {
    //             if ((ol_string_graph->ver[i].in_edge[j].adjvex&1) == 0) fprintf(fp_g, "%d.B,", ol_string_graph->ver[i].in_edge[j].adjvex>>1);
    //             else if ((ol_string_graph->ver[i].in_edge[j].adjvex&1) == 1)    fprintf(fp_g, "%d.E,", ol_string_graph->ver[i].in_edge[j].adjvex>>1);
    //             if ((i&1) == 0) fprintf(fp_g, "%d.B,", i>>1);
    //             else if ((i&1) == 1)    fprintf(fp_g, "%d.E,", i>>1);
    //             fprintf(fp_g, "%d,", ol_string_graph->ver[i].in_edge[j].end_pos);
    //             if (ol_string_graph->ver[i].in_edge[j].str == 0)
    //                 fprintf(fp_g, "%d,%d,+%d,", ol_string_graph->ver[i].in_edge[j].reduce, ol_string_graph->ver[i].in_edge[j].id, ol_string_graph->ver[i].in_edge[j].utg_id);
    //             else if (ol_string_graph->ver[i].in_edge[j].str == 1)
    //                 fprintf(fp_g, "%d,%d,-%d,", ol_string_graph->ver[i].in_edge[j].reduce, ol_string_graph->ver[i].in_edge[j].id, ol_string_graph->ver[i].in_edge[j].utg_id);
    //             fprintf(fp_g, "%d,", ol_string_graph->ver[i].in_edge[j].edge_list.ver_n);
    //             for (size_t k = 0; k < ol_string_graph->ver[i].in_edge[j].edge_list.ver_n; k++)
    //             {
    //                 fprintf(fp_g, "[%d %d %d-%d,%d-%d] ", ol_string_graph->ver[i].in_edge[j].edge_list.vex[k]&1, ol_string_graph->ver[i].in_edge[j].edge_list.vex[k] >> 1, ol_string_graph->ver[i].in_edge[j].edge_list.edge_sta[k], ol_string_graph->ver[i].in_edge[j].edge_list.edge_end[k], ol_string_graph->ver[i].in_edge[j].edge_list.read_sta[k], ol_string_graph->ver[i].in_edge[j].edge_list.read_end[k]);
    //             }
    //             fprintf(fp_g, "\n");
    //         }
    //     }
    // }

    fclose(fp_g);
#endif

    if (vi != NULL) free(vi);

    return n_utg;
}