#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <assert.h>

#include "bseq.h"
#include "paf.h"
#include "overlapping.h"
#include "layout.h"
// #include "linking.h"
// #include "assemblygraph.h"

uint32_t load_oves_info(const char *temp_ove_dir, OVE_C_t *ove_cl)
{
    paf_file_t *fp;
    paf_rec_t r;
    uint32_t vertex_num = 0;

    fp = paf_open(temp_ove_dir);
    if (!fp)
    {
        fprintf(stderr, "[%s] could not open PAF file %s\n", __func__, temp_ove_dir);
        exit(1);
    }

    while (paf_read(fp, &r) >= 0) //r: a paf row
    {
        if (ove_cl->n == ove_cl->m)
        {
            ove_cl->m = ove_cl->m << 1;
            ove_cl->ove = (OVE_t *)realloc(ove_cl->ove, ove_cl->m * sizeof(OVE_t));
        }
        ove_cl->ove[ove_cl->n].qid = r.qid; ove_cl->ove[ove_cl->n].tid = r.tid;
        ove_cl->ove[ove_cl->n].ql = r.ql; ove_cl->ove[ove_cl->n].tl = r.tl;
        ove_cl->ove[ove_cl->n].qs = r.qs; ove_cl->ove[ove_cl->n].qe = r.qe;
        ove_cl->ove[ove_cl->n].ts = r.ts; ove_cl->ove[ove_cl->n].te = r.te;
        ove_cl->ove[ove_cl->n].rev = r.rev; ove_cl->ove[ove_cl->n].mbp = r.mb; ove_cl->ove[ove_cl->n++].mln = r.ml;
        if (r.qid > vertex_num){
            vertex_num = r.qid;
            // printf("%ld:%d\t",ove_cl->n,vertex_num);
        }
    }
    // printf("\n");

    paf_close(fp);
    return vertex_num;
}

int layout_workflow(const char *read_fastq, const char *temp_ove_dir, const char *temp_ove_link_dir)
{
    uint32_t vertex_num;
    OVE_C_t *ove_cl;
    ove_cl = (OVE_C_t *)calloc(1, sizeof(OVE_C_t));
    ove_cl->n = 0;
    ove_cl->m = 65535;
    ove_cl->ove = (OVE_t*)calloc(ove_cl->m, sizeof(OVE_t));
    if (ove_cl->ove == NULL)
    {
        fprintf(stderr, "[%s] could alloc ove_cl->ove\n", __func__);
        exit(1);
    }

    OVE_C_t *ove_link;
    ove_link = (OVE_C_t *)calloc(1, sizeof(OVE_C_t));
    ove_link->n = 0;
    ove_link->m = 655350;
    ove_link->ove = (OVE_t*)calloc(ove_link->m, sizeof(OVE_t));
    if (ove_link->ove == NULL)
    {
        fprintf(stderr, "[%s] could alloc ove_link->ove\n", __func__);
        exit(1);
    }

    vertex_num = load_oves_info(temp_ove_dir, ove_cl);
    load_oves_info(temp_ove_link_dir, ove_link);

    fprintf(stderr, "[Layout] Loading paf files, %ld index link, %ld query link, inital vertex num %d...\n", ove_cl->n, ove_link->n, vertex_num);

    StringGraph *ol_string_graph = NULL;
    ol_string_graph = (StringGraph *)calloc(1, sizeof(StringGraph));
    if (ol_string_graph == NULL)
    {
        fprintf(stderr, "[%s] Failed to allocate graph memory!\n", __func__);
        exit(1);
    }



    if (ol_string_graph != NULL)   free(ol_string_graph);

    if (ove_cl->ove != NULL) free(ove_cl->ove);
    if (ove_link->ove != NULL) free(ove_link->ove);

    return 0;
}