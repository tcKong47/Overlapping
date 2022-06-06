#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#include "bseq.h"
#include "ktime.h"
#include "index.h"
#include "ksort.h"
#include "thread.h"
#include "bit_operation.h"
#include "overlapping.h"

#define sort_key_hash(h) ((h).x)
KRADIX_SORT_INIT(mi, MINIMIZER_t, sort_key_hash, sizeof(uint64_t))

uint32_t sketching_core(READ_t *seq, int ridx, MINIMIZER_t *seq_mi, int w, int k)
{
    uint32_t i, j, buf_pos = 0, min_pos = 0, len = 0, seq_mi_n = 0;
    uint64_t kmer[2] = {0, 0}, shift, mask;
    MINIMIZER_t temp_mi = {UINT64_MAX, UINT64_MAX}, min_mi = {UINT64_MAX, UINT64_MAX};
    MINIMIZER_t *buf = NULL;

    shift = 2 * (k - 1);        // 2 * 14
    mask = (1ULL << 2 * k) - 1; // (111...111) 30
    buf = (MINIMIZER_t *)calloc(w, sizeof(MINIMIZER_t));
	if (buf == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB buf memory error!\n", __func__, w * sizeof(MINIMIZER_t) / 1024 /1024 / 1024);
		exit(1);
	}

    for (i = 0; i < seq->read_length; ++i)
    {
        temp_mi.x = UINT64_MAX;
        temp_mi.y = UINT64_MAX;
        int c = seq_nt4_table[(uint8_t)seq->read_seq[i]];
        if (c < 4)
        {
            int z;
            kmer[0] = (kmer[0] << 2 | c) & mask;            // forward k-mer
            kmer[1] = (kmer[1] >> 2) | (3ULL - c) << shift; // reverse k-mer
            if (kmer[0] == kmer[1]) continue;
            z = kmer[0] < kmer[1] ? 0 : 1; // strand
            if (++len >= k)
            {
                temp_mi.x = hash64(kmer[z], mask);
				temp_mi.y = (uint64_t)ridx << 32 | (uint32_t)i << 1 | z;
            }
        }
        else //if meet char 'N', skip this region
        {
            len = 0;
        }

        buf[buf_pos] = temp_mi;

        if (len == w + k - 1) // for the first window
        {
            seq_mi[seq_mi_n++] = min_mi;
        }

        if (temp_mi.x <= min_mi.x)
        {
            min_mi = temp_mi;
            min_pos = buf_pos;
            if (len >= w + k)
            {
				seq_mi[seq_mi_n++] = min_mi;
			}
        }
        else if (min_pos == buf_pos) //min_mi has been removed outside the window
        {
            min_mi.x = UINT64_MAX;
            for (j = 0; j < w; ++j)
            {
                if (buf[buf_pos].x < min_mi.x)
                {
                    min_mi.x = buf[buf_pos].x;
                    min_mi.y = buf[buf_pos].y;
					seq_mi[seq_mi_n++] = min_mi;
				}
            }
        }

        if (++buf_pos == w) buf_pos = 0;
    }
    // printf("\n");

    if (buf != NULL) {free(buf); buf = NULL;}

    return seq_mi_n;
}

uint32_t filter_rep_mi(mini_idx_t *mi)
{
	uint32_t i, rep_idx, single_n = 0;
	uint32_t *rep_kmer, rep_tmp = 1, rep_n = 0;
	if (r_n <= 0.) return UINT32_MAX;

	rep_kmer = (uint32_t *)calloc(mi->mm.mi_n, sizeof(uint32_t));
	if (rep_kmer == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB rep_kmer memory error!\n", __func__, mi->mm.mi_n * sizeof(uint32_t) / 1024 /1024 / 1024);
		exit(1);
	}

	uint64_t mi_x_tmp = mi->mm.mi_v[0].x;
	for (i = 1; i < mi->mm.mi_n; i++)
	{
		if (mi->mm.mi_v[i].x == mi_x_tmp) rep_tmp++;
		else
		{
			if (rep_tmp > 1) rep_kmer[rep_n++] = rep_tmp;
			else if (rep_tmp == 1) single_n += 1;
			rep_tmp = 1;
			mi_x_tmp = mi->mm.mi_v[i].x;
		}
	}
	if (rep_tmp > 1) rep_kmer[rep_n++] = rep_tmp; // last one

	if (rep_n == 0)
	{
		rep_tmp = 1;
		fprintf(stderr, "[Index] %d distinct k-mer..., ignore k-mer which appearing times > %d...\n", single_n, rep_tmp);
	}
	else
	{
		qsort(rep_kmer, rep_n, sizeof(uint32_t), compare_rep);
		rep_idx = rep_n * r_n;
		rep_tmp = rep_kmer[rep_idx];
		fprintf(stderr, "[Index] %d distinct k-mer..., ignore k-mer which appearing times > %d...\n", rep_n + single_n, rep_tmp);
	}

	if (rep_kmer != NULL) {free(rep_kmer); rep_kmer = NULL;}

	return rep_tmp;
}

// shared data of all threads
typedef struct
{
	uint32_t index_batch_size;
	uint32_t query_batch_size;
	bseq_file_t *fp;
	uint64_t n_base;
	uint32_t n_processed;
	uint32_t n_sta;
	mini_idx_t *mi;
	uint32_t *mi_num;
	uint32_t n_sta_idx;
	uint32_t *iter_idx;
} pipeline_t;

typedef struct
{
	uint32_t n_sta_idx;
	uint32_t n_seq;
	READ_t *read_info;
} step_index_t;

static void *index_pipeline(void *shared, int step, void *in)
{
	int32_t i, j;
	pipeline_t *p = (pipeline_t *)shared;
	if (step == 0)
	{
		step_index_t *s;
		if (p->n_base > p->query_batch_size) return 0;
		s = (step_index_t *)calloc(1, sizeof(step_index_t));
		if (s == NULL)
		{
			fprintf(stderr, "[%s] calloc s memory error!\n", __func__);
			exit(1);
		}

		s->n_sta_idx = p->n_sta_idx;
		s->read_info = bseq_read_idx(p->fp, p->index_batch_size, &s->n_seq, &p->n_sta, &p->n_sta_idx, p->iter_idx);

		if (s->n_seq)
		{
			p->n_processed += s->n_seq;
			for (i = 0; i < s->n_seq; i++)
			{
				p->n_base += s->read_info[i].read_length;
			}
			return s;
		}
		else
		{
			if (s != NULL) {free(s); s = NULL;}
		}
	}
	else if (step == 1)
	{
		step_index_t *s = (step_index_t *)in;
		uint32_t seq_mi_n = 0, key;
		MINIMIZER_t *seq_mi = NULL;
		for (i = 0; i < s->n_seq; i++)
		{
			READ_t *read = &s->read_info[i];
			seq_mi = (MINIMIZER_t *)realloc(seq_mi, read->read_length * sizeof(MINIMIZER_t));
			if (seq_mi == NULL)
			{
				fprintf(stderr, "[%s] calloc %ldGB seq_mi memory error!\n", __func__, read->read_length * sizeof(MINIMIZER_t) / 1024 /1024 / 1024);
				exit(1);
			}
			seq_mi_n = sketching_core(&s->read_info[i], i+s->n_sta_idx, seq_mi, p->mi->w, p->mi->k);
			for (j = 0; j < seq_mi_n; j++)
			{
				if (p->mi->mm.mi_n == p->mi->mm.mi_m)
				{
					p->mi->mm.mi_m = p->mi->mm.mi_m << 1;
					p->mi->mm.mi_v = (MINIMIZER_t *)realloc(p->mi->mm.mi_v, p->mi->mm.mi_m * sizeof(MINIMIZER_t));
					if (p->mi->mm.mi_v == NULL)
					{
						fprintf(stderr, "[%s] calloc %ldGB p->mi->mm.mi_v memory error!\n", __func__, p->mi->mm.mi_m * sizeof(MINIMIZER_t) / 1024 /1024 / 1024);
						exit(1);
					}
				}
				p->mi->mm.mi_v[p->mi->mm.mi_n++] = seq_mi[j];
				if (p->mi->mm.mi_n > p->mi->mm.mi_m) printf("[%s]seq_mi memory leak..,%ld > %ld\n",__func__,p->mi->mm.mi_n,p->mi->mm.mi_m);
				key = (seq_mi[j].x >> (k - l) * 2) & p->mi->mi_mask;
				p->mi_num[key]++;
			}
			if (read->read_seq != NULL) {free(read->read_seq); read->read_seq = NULL;}
			if (read->read_name != NULL) {free(read->read_name); read->read_name = NULL;}
		}
		if (seq_mi != NULL) {free(seq_mi); seq_mi = NULL;}
		if (s->read_info != NULL) {free(s->read_info); s->read_info = NULL;}
		if (s != NULL) {free(s); s = NULL;}
	}

	return 0;
}

mini_idx_t *indexing_input_read(bseq_file_t *fp, uint32_t *rep_n, uint32_t *n_sta, uint32_t *in_idx_n_sta, read_stat_t *read_stat)
{
	pipeline_t pl;

	pl.index_batch_size = batch_size_base / 5;
	pl.query_batch_size = batch_size_base / 1; // 33000+ reads, 24000+bp
	pl.fp = fp;
	pl.n_base = 0;
	pl.n_processed = 0;
	pl.n_sta = *n_sta;
	pl.n_sta_idx = *in_idx_n_sta;
	pl.iter_idx = read_stat->iter_idx;

	pl.mi = (mini_idx_t *)calloc(1, sizeof(mini_idx_t));
	pl.mi->k = k;
	pl.mi->l = l;
	pl.mi->w = w;
	pl.mi->bucket_num = 1 << (pl.mi->l << 1);	   // bucket num
	pl.mi->mi_mask = (1ULL << 2 * pl.mi->l) - 1; // (111...111) 20
	pl.mi->mm.mi_n = 0;
	pl.mi->mm.mi_m = pl.mi->bucket_num;
	pl.mi->mm.mi_v = (MINIMIZER_t *)calloc(pl.mi->mm.mi_m, sizeof(MINIMIZER_t));
	if (pl.mi->mm.mi_v == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB pl.mi->mm.mi_v memory error!\n", __func__, pl.mi->mm.mi_m * sizeof(MINIMIZER_t) / 1024 /1024 / 1024);
		exit(1);
	}
	pl.mi->mi_count = (uint32_t *)calloc(pl.mi->bucket_num + 1, sizeof(uint32_t));
	if (pl.mi->mi_count == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB pl.mi->mi_count memory error!\n", __func__, pl.mi->bucket_num * sizeof(uint32_t) / 1024 /1024 / 1024);
		exit(1);
	}
	pl.mi_num = (uint32_t *)calloc(pl.mi->bucket_num, sizeof(uint32_t));
	if (pl.mi_num == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB pl.mi_num memory error!\n", __func__, pl.mi->bucket_num * sizeof(uint32_t) / 1024 /1024 / 1024);
		exit(1);
	}

	kt_pipeline(thread_n, index_pipeline, &pl, 2);
	if (pl.mi->mm.mi_n == 0)
	{
		if (pl.mi_num != NULL) {free(pl.mi_num); pl.mi_num = NULL;}
		return pl.mi;
	}
	fprintf(stderr, "[Index] Total %d indexed reads, %ld minimizers(%ld); step: %f; total memory %.3f GB; soring minimizers...\n", pl.n_processed, pl.mi->mm.mi_n, pl.mi->mm.mi_m, (double)pl.n_base / pl.mi->mm.mi_n, pl.mi->mm.mi_n * sizeof(MINIMIZER_t) / 1024.0 / 1024.0 / 1024.0);

	pl.mi->mi_count[0] = 0;
	for (int i = 1; i < pl.mi->bucket_num; ++i)
	{
		pl.mi->mi_count[i] = pl.mi->mi_count[i - 1] + pl.mi_num[i - 1];
	}
	pl.mi->mi_count[pl.mi->bucket_num] = pl.mi->mm.mi_n;

	radix_sort_mi(&(pl.mi->mm.mi_v[0]), &(pl.mi->mm.mi_v[pl.mi->mm.mi_n - 1]) + 1);

	*rep_n = filter_rep_mi(pl.mi);
	*n_sta = pl.n_sta;
	*in_idx_n_sta = pl.n_sta_idx;
	fprintf(stderr, "[%s]\tReal time:%.3f sec\tCPU:%.3f sec\tMemory peak:%.3f GB\n", __func__, realtime() - realtime0, cputime(), peak_memory() / 1024.0 / 1024.0);

	if (pl.mi_num != NULL) {free(pl.mi_num); pl.mi_num = NULL;}

	return pl.mi;
}