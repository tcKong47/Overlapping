#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <assert.h>
#include <math.h>

#include "bseq.h"
#include "kvec.h"
#include "ktime.h"
#include "index.h"
#include "graph.h"
#include "thread.h"
#include "overlapping.h"
#include "bit_operation.h"

// shared data of all threads
typedef struct
{
	uint32_t iter;
	read_stat_t *read_stat;
	uint32_t max_len;
	uint32_t x_len;
	uint32_t rep_n;
	uint32_t batch_i;
	uint32_t n_threads;
	uint32_t n_processed;
	uint32_t n_total_read;
	uint32_t index_batch_size;
	uint32_t query_batch_size;
	bseq_file_t *fp;
	mini_idx_t *mi;
	uint32_t top_n;
	uint32_t index_n_sta;     // for iter_idx, start rid of current batch
	uint32_t index_idx_n_sta; // for iter_idx, start index of current batch
	uint32_t index_idx_n_end; // for iter_idx, end index of current batch, start index of next batch
	uint32_t map_n_sta;     // for iter_map, start rid of current batch
	uint32_t map_idx_n_sta; // for iter_map, start index of current batch
	uint32_t map_idx_n_end; // for iter_map, end index of current batch, start index of next batch
	OVE_C_t *ove_cl;
	int m_b;
	int min_ove;
	int max_hang;
	float max_hang_ratio;
} pipeline_t;

typedef struct
{
	uint32_t n_seq;
	uint32_t *len;
} step_sort_t;

typedef struct
{
	const pipeline_t *p;
	uint32_t n_seq;
	uint32_t n_sta;
	uint32_t n_sta_idx;
	READ_t *read_info;
	thread_buf_t **buf;
	uint32_t ove_alloc;
} step_query_t;

double mapping_time = 0., cpu_time = 0.;

static void *sort_read_length_pipeline(void *shared, int step, void *in)
{
	int i;
	pipeline_t *p = (pipeline_t *)shared;
	if (step == 0)
	{
		step_sort_t *s = NULL;
		s = (step_sort_t *)calloc(1, sizeof(step_sort_t));
		if (s == NULL)
		{
			fprintf(stderr, "[%s] calloc s memory error!\n", __func__);
			exit(1);
		}
		s->len = bseq_read_l(p->fp, p->index_batch_size, &s->n_seq);
		if (s->n_seq)
		{
			p->read_stat->rlen = (uint32_t *)realloc(p->read_stat->rlen, (p->n_processed + s->n_seq) * sizeof(uint32_t));
			if (p->read_stat->rlen == NULL)
			{
				fprintf(stderr, "[%s] calloc %ldGB p->read_stat->rlen memory error!\n", __func__, (p->n_processed + s->n_seq) * sizeof(uint32_t) / 1024 / 1024 / 1024);
				exit(1);
			}
			for (i = 0; i < s->n_seq; i++)
			{
				p->read_stat->rlen[p->n_processed++] = s->len[i];
				if (s->len[i] > p->max_len) p->max_len = s->len[i];
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
		step_sort_t *s = (step_sort_t *)in;
		qsort(s->len, s->n_seq, sizeof(uint32_t), compare_rlen);
		uint32_t read_n = s->n_seq * x_read * 0.01;
		p->x_len += s->len[read_n - 1];
		p->batch_i += 1;
		if (s->len != NULL) {free(s->len); s->len = NULL;}
		if (s != NULL) {free(s); s = NULL;}
	}
	return 0;
}

int binary_search(MINIMIZER_t minmer, int k, int l, uint64_t *pos, mini_idx_t *mi)
{
	uint64_t key;
	uint64_t pos_s = 0, pos_e = 0, pos_m = 0;
	uint64_t mask_l = (1ULL << 2 * l) - 1; // (111...111) 20

	pos[0] = pos[1] = -1;
	key = (minmer.x >> (k-l)*2) & mask_l;
	pos_s = mi->mi_count[key];
	pos_e = mi->mi_count[key+1] - 1;

	if (minmer.x < mi->mm.mi_v[pos_s].x || minmer.x > mi->mm.mi_v[pos_e].x)	return -1;
	while(pos_s <= pos_e)
	{
		pos_m = (pos_s + pos_e)/2;
		if (mi->mm.mi_v[pos_m].x == minmer.x)
		{
			pos[0] = pos[1] = pos_m;
			//left extend
			int sub_s = pos_s, sub_e = pos_m - 1, sub_m;
			while(sub_s <= sub_e)
			{
				sub_m = (sub_s + sub_e) / 2;
				if (mi->mm.mi_v[sub_m].x == minmer.x)
				{
					pos[0] = sub_m;
					sub_e = sub_m - 1;
				}
				else if (mi->mm.mi_v[sub_m].x < minmer.x)	sub_s = sub_m + 1;
				else if (mi->mm.mi_v[sub_m].x > minmer.x)	exit(1);
			}
			//right extend
			sub_s = pos_m + 1;
			sub_e = pos_e;
			while(sub_s <= sub_e)
			{
				sub_m = (sub_s + sub_e) / 2;
				if (mi->mm.mi_v[sub_m].x == minmer.x)
				{
					pos[1] = sub_m;
					sub_s = sub_m + 1;
				}
				else if (mi->mm.mi_v[sub_m].x > minmer.x)	sub_e = sub_m -1 ;
				else if (mi->mm.mi_v[sub_m].x < minmer.x)	exit(1);
			}
			return 1;
		}
		else if (mi->mm.mi_v[pos_m].x > minmer.x)	pos_e = pos_m - 1;
		else if (mi->mm.mi_v[pos_m].x < minmer.x)	pos_s = pos_m + 1;
	}
	return -1;
}

uint32_t finding_hits(void *shared_data, int64_t qidx, uint32_t pid)
{
	int result;
	step_query_t *data = (step_query_t *)shared_data;
	READ_t *query = &data->read_info[qidx];
	const pipeline_t *p = data->p;
	mini_idx_t *mi = data->p->mi;
	mini_t *seq_mi = data->buf[pid]->seq_mi;
	uint32_t v_n = 0, *v_m = &data->buf[pid]->v_m;
	uint32_t i, ii, hits_num = 0;
	uint64_t range[2] = {UINT64_MAX, UINT64_MAX};
	uint32_t pos_q, pos_t, strand_mask = 1ULL, pos_mask = (1ULL << 32) - 1;
	uint32_t su_i = 0, j, tid_tmp;
	uint32_t s1, e1, Eindel = 20, tstr_tmp, qstr_tmp, cov_tmp;
	uint32_t rid;

	assert(query->read_length <= seq_mi->mi_m);
	seq_mi->mi_n = sketching_core(&data->read_info[qidx], qidx, seq_mi->mi_v, mi->w, mi->k);
	if (seq_mi->mi_n >= seq_mi->mi_m) printf("[%s] use unlegal memory, %ld-%ld, something wrong with it...", __func__, seq_mi->mi_n, seq_mi->mi_m);

	for (i = 0; i < seq_mi->mi_n; ++i)
	{
		result = binary_search(seq_mi->mi_v[i], p->mi->k, p->mi->l, range, mi);
		if (result != -1)
		{
			hits_num = range[1] - range[0] + 1;
			if (hits_num > 0 && hits_num <= p->rep_n) //meeting repeat region, drop these hits
			{
				for (ii = range[0]; ii <= range[1]; ++ii)
				{
					rid = p->read_stat->iter_idx[mi->mm.mi_v[ii].y >> 32];
					if (rid != query->rid)
					{
						pos_q = ((seq_mi->mi_v[i].y & pos_mask) >> 1);
						pos_t = ((mi->mm.mi_v[ii].y & pos_mask) >> 1);

						if (v_n >= (*v_m))
						{
							(*v_m) = (*v_m) == 0 ? 256 : (*v_m) << 1;
							data->buf[pid]->vertex_mr = (MR_t *)realloc(data->buf[pid]->vertex_mr, (*v_m) * sizeof(MR_t));
							if (data->buf[pid]->vertex_mr == NULL)
							{
								fprintf(stderr, "[%s] calloc %ldGB data->buf[pid]->vertex_mr memory error!\n", __func__, (*v_m) * sizeof(MR_t) / 1024 / 1024 / 1024);
								exit(1);
							}
							printf("[%s] will not use these code, something wrong with it...",__func__);
						}
						data->buf[pid]->vertex_mr[v_n].qs = pos_q + 1 - k;
						data->buf[pid]->vertex_mr[v_n].qe = pos_q;
						data->buf[pid]->vertex_mr[v_n].ts = pos_t + 1 - k;
						data->buf[pid]->vertex_mr[v_n].te = pos_t;
						data->buf[pid]->vertex_mr[v_n].cov = k;
						data->buf[pid]->vertex_mr[v_n].t_id = mi->mm.mi_v[ii].y >> 32;
						data->buf[pid]->vertex_mr[v_n].tstr = (mi->mm.mi_v[ii].y & strand_mask);
						data->buf[pid]->vertex_mr[v_n++].qstr = (seq_mi->mi_v[i].y & strand_mask);
					}
				}
			}
		}
	}

	data->buf[pid]->v_n = v_n;
	if (v_n == 0)	return v_n;

	// merge co-linear match region in the same target read
	qsort(data->buf[pid]->vertex_mr, v_n, sizeof(MR_t), compare_tid);

	tid_tmp = data->buf[pid]->vertex_mr[0].t_id;
	tstr_tmp = data->buf[pid]->vertex_mr[0].tstr;
	qstr_tmp = data->buf[pid]->vertex_mr[0].qstr;
	j = 0;
	while (j < v_n)
	{
		s1 = j++;
		cov_tmp = data->buf[pid]->vertex_mr[s1].cov;
		while ((tid_tmp == data->buf[pid]->vertex_mr[j].t_id) && ((tstr_tmp ^ data->buf[pid]->vertex_mr[j].tstr) == (qstr_tmp ^ data->buf[pid]->vertex_mr[j].qstr)) && (data->buf[pid]->vertex_mr[j].qs >= data->buf[pid]->vertex_mr[j - 1].qs) && (j < v_n))
		{
			int diff = (int)(data->buf[pid]->vertex_mr[j].qs - data->buf[pid]->vertex_mr[j - 1].qe - 1);
			int diff_q = (int)(data->buf[pid]->vertex_mr[j].qs - data->buf[pid]->vertex_mr[j - 1].qs);
			int diff_t = data->buf[pid]->vertex_mr[j].tstr == data->buf[pid]->vertex_mr[j].qstr ? (int)(data->buf[pid]->vertex_mr[j].ts - data->buf[pid]->vertex_mr[j - 1].ts) : (int)(data->buf[pid]->vertex_mr[j - 1].te - data->buf[pid]->vertex_mr[j].te);
			// printf("s1 %d,%d %d, diff %d %d, %d, %d, %d\n",s1,j-1,j,diff,waitingLen,diff_q,diff_t,abs(diff_t - diff_q));
			if (diff > waitingLen)
				break;
			if (abs(diff_t - diff_q) < Eindel)
			{
				cov_tmp += (diff > 0) ? data->buf[pid]->vertex_mr[j++].cov : (diff + data->buf[pid]->vertex_mr[j++].cov);
			}
			else
				break;
		}
		e1 = j - 1;

		data->buf[pid]->vertex_mr[su_i].t_id = data->buf[pid]->vertex_mr[s1].t_id;
		data->buf[pid]->vertex_mr[su_i].qs = data->buf[pid]->vertex_mr[s1].qs;
		data->buf[pid]->vertex_mr[su_i].qe = data->buf[pid]->vertex_mr[e1].qe;
		data->buf[pid]->vertex_mr[su_i].ts = data->buf[pid]->vertex_mr[s1].tstr == data->buf[pid]->vertex_mr[s1].qstr ? data->buf[pid]->vertex_mr[s1].ts : data->buf[pid]->vertex_mr[e1].ts;
		data->buf[pid]->vertex_mr[su_i].te = data->buf[pid]->vertex_mr[e1].tstr == data->buf[pid]->vertex_mr[e1].qstr ? data->buf[pid]->vertex_mr[e1].te : data->buf[pid]->vertex_mr[s1].te;
		data->buf[pid]->vertex_mr[su_i].tstr = data->buf[pid]->vertex_mr[e1].tstr;
		data->buf[pid]->vertex_mr[su_i].qstr = data->buf[pid]->vertex_mr[e1].qstr;
		data->buf[pid]->vertex_mr[su_i++].cov = cov_tmp;

		tid_tmp = data->buf[pid]->vertex_mr[j].t_id;
		tstr_tmp = data->buf[pid]->vertex_mr[j].tstr;
		qstr_tmp = data->buf[pid]->vertex_mr[j].qstr;
		if (su_i >= data->buf[pid]->v_m) printf("[%s] use unlegal memory, %d-%d, something wrong with it...", __func__, su_i, data->buf[pid]->v_m);
	}
	data->buf[pid]->v_n = su_i;

	uint32_t count = 1, su_ii = 0;
	cov_tmp = data->buf[pid]->vertex_mr[0].cov;
	tid_tmp = data->buf[pid]->vertex_mr[0].t_id;
	for (j = 1; j < su_i; ++j)
	{
		if (data->buf[pid]->vertex_mr[j].t_id != tid_tmp)
		{
			if (count == 1 && cov_tmp <= k)
			{
				count = 1;
				cov_tmp = data->buf[pid]->vertex_mr[j].cov;
				tid_tmp = data->buf[pid]->vertex_mr[j].t_id;
			}
			else
			{
				// printf("%d %d %d %d\n",count,cov_tmp,tid_tmp,j);
				for (int k = count; k > 0; k--)
				{
					// printf("%d %d %d %d\n", su_ii, j, k, j-k);
					data->buf[pid]->vertex_mr[su_ii].t_id = data->buf[pid]->vertex_mr[j - k].t_id;
					data->buf[pid]->vertex_mr[su_ii].qs = data->buf[pid]->vertex_mr[j - k].qs;
					data->buf[pid]->vertex_mr[su_ii].qe = data->buf[pid]->vertex_mr[j - k].qe;
					data->buf[pid]->vertex_mr[su_ii].ts = data->buf[pid]->vertex_mr[j - k].ts;
					data->buf[pid]->vertex_mr[su_ii].te = data->buf[pid]->vertex_mr[j - k].te;
					data->buf[pid]->vertex_mr[su_ii].tstr = data->buf[pid]->vertex_mr[j - k].tstr;
					data->buf[pid]->vertex_mr[su_ii].qstr = data->buf[pid]->vertex_mr[j - k].qstr;
					data->buf[pid]->vertex_mr[su_ii++].cov = data->buf[pid]->vertex_mr[j - k].cov;
				}
				count = 1;
				cov_tmp = data->buf[pid]->vertex_mr[j].cov;
				tid_tmp = data->buf[pid]->vertex_mr[j].t_id;
			}
		}
		else
		{
			cov_tmp += data->buf[pid]->vertex_mr[j].cov;
			count++;
		}
		if (su_ii >= data->buf[pid]->v_m) printf("[%s] use unlegal memory, %d-%d, something wrong with it...", __func__, su_ii, data->buf[pid]->v_m);
	}

	// the last one
	if (!(count == 1 && cov_tmp <= k))
	{
		for (int k = count; k > 0; k--)
		{
			data->buf[pid]->vertex_mr[su_ii].t_id = data->buf[pid]->vertex_mr[j - k].t_id;
			data->buf[pid]->vertex_mr[su_ii].qs = data->buf[pid]->vertex_mr[j - k].qs;
			data->buf[pid]->vertex_mr[su_ii].qe = data->buf[pid]->vertex_mr[j - k].qe;
			data->buf[pid]->vertex_mr[su_ii].ts = data->buf[pid]->vertex_mr[j - k].ts;
			data->buf[pid]->vertex_mr[su_ii].te = data->buf[pid]->vertex_mr[j - k].te;
			data->buf[pid]->vertex_mr[su_ii].tstr = data->buf[pid]->vertex_mr[j - k].tstr;
			data->buf[pid]->vertex_mr[su_ii].qstr = data->buf[pid]->vertex_mr[j - k].qstr;
			data->buf[pid]->vertex_mr[su_ii++].cov = data->buf[pid]->vertex_mr[j - k].cov;
		}
	}

	if (su_ii >= data->buf[pid]->v_m) printf("[%s] use unlegal memory, %d-%d, something wrong with it...", __func__, su_ii, data->buf[pid]->v_m);
	data->buf[pid]->v_n = su_ii;

	// printf("qid\thid\tcov\tqstrand\tqs\tqe\ttid\ttstrand\tts\tte\tlen\n");
	// for (int vi = 0; vi < su_ii; ++vi)
	// {
	// 	printf("%d\t%d\t%d\t%d\t%ld\t%ld\t%d\t%d\t%ld\t%ld\n", qidx, vi, data->buf[pid]->vertex_mr[vi].cov, data->buf[pid]->vertex_mr[vi].qstr, data->buf[pid]->vertex_mr[vi].qs, data->buf[pid]->vertex_mr[vi].qe, data->buf[pid]->vertex_mr[vi].t_id, data->buf[pid]->vertex_mr[vi].tstr, data->buf[pid]->vertex_mr[vi].ts, data->buf[pid]->vertex_mr[vi].te);
	// }
	// printf("\n");

	return su_ii;
}

void cp_ove(OVE_t *dest_ove, const OVE_t *src_ove)
{
	dest_ove->qid = src_ove->qid;
	dest_ove->ql = src_ove->ql;
	dest_ove->qs = src_ove->qs;
	dest_ove->qe = src_ove->qe;
	dest_ove->tid = src_ove->tid;
	dest_ove->tl = src_ove->tl;
	dest_ove->ts = src_ove->ts;
	dest_ove->te = src_ove->te;
	dest_ove->rev = src_ove->rev;
	dest_ove->mbp = src_ove->mbp;
	dest_ove->mln = src_ove->mln;
	dest_ove->score = src_ove->score;
}

int get_overlap_info(step_query_t *data, thread_buf_t *buf, uint32_t q_idx, uint32_t pid)
{
	READ_t *query = &data->read_info[q_idx];
	int32_t j;
	uint32_t i = 0, matching_bases = 0;
	OVE_t ove_tmp;
	OVE_t *ove_cl_tmp = NULL;
	uint32_t n_tmp = 0;
	uint32_t qpre = 0, qsuf = 0, tpre = 0, tsuf = 0;
	uint32_t left_overhang = 0, right_overhang = 0, overhang_tr = 0, overhang = 0;
	uint32_t max_index_n = buf->max_index_n;
	uint32_t *max_index = buf->max_index;
	PATH_t *dist_path = buf->path;
	MR_t *vertex_mr = buf->vertex_mr;

	if (max_index_n == 0) return 0;
	ove_cl_tmp = (OVE_t*)calloc(max_index_n,sizeof(OVE_t));
	if (ove_cl_tmp == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB ove_cl_tmp memory error!\n", __func__, max_index_n * sizeof(OVE_t) / 1024 / 1024 / 1024);
		exit(1);
	}
	for (i = 0; i < max_index_n; ++i)
	{
		j = (int32_t)max_index[i];
		matching_bases = 0;
		ove_tmp.score = dist_path[j].dist;
		// first iter vertex_mr.t_id is the rid of index read; second iter vertex_mr.t_id is the index of index read in array read_stat->iter_idx
		ove_tmp.tid  = vertex_mr[j].t_id;
		ove_tmp.qid = query->rid;
		ove_tmp.ql = query->read_length;
		ove_tmp.tl = data->p->read_stat->rlen[data->p->read_stat->iter_idx[ove_tmp.tid]];
		ove_tmp.rev = vertex_mr[j].tstr == vertex_mr[j].qstr ? 0 : 1;
		ove_tmp.qe = vertex_mr[j].qe;
		ove_tmp.te = ove_tmp.rev == 0 ? vertex_mr[j].te : vertex_mr[j].ts;

		while(j != -1)
		{
			ove_tmp.qs = vertex_mr[j].qs;
			ove_tmp.ts = ove_tmp.rev == 0 ? vertex_mr[j].ts : vertex_mr[j].te;
			matching_bases += vertex_mr[j].cov;
			j = (int32_t)dist_path[j].pre_node;
			if (j >= (int32_t)buf->v_n) printf("[%s] 1 use unlegal memory, %d-%d, something wrong with it...", __func__, j, buf->v_n);
		}

		if (ove_tmp.rev == 0)
		{
			qpre = ove_tmp.qs; qsuf = ove_tmp.ql - ove_tmp.qe;
			tpre = ove_tmp.ts; tsuf = ove_tmp.tl - ove_tmp.te;
		}
		else if (ove_tmp.rev == 1)
		{
			qpre = ove_tmp.qs; qsuf = ove_tmp.ql - ove_tmp.qe;
			tsuf = ove_tmp.te; tpre = ove_tmp.tl - ove_tmp.ts;
		}

		left_overhang  = qpre < tpre ? qpre : tpre;
		right_overhang = qsuf < tsuf ? qsuf : tsuf;
		overhang = left_overhang + right_overhang;
		ove_tmp.mbp = matching_bases;
		ove_tmp.mln = abs(ove_tmp.te - ove_tmp.ts) + 1;

		// if (query->rid == 37) printf("read %d,%d: %d\t%d\t%d\t%d\t%.6f\n", ove_tmp.qid, data->p->read_stat->iter_idx[ove_tmp.tid], ove_tmp.mbp, ove_tmp.mln, left_overhang, right_overhang, ove_tmp.score);
		// overhang_tr = data->p->max_hang;
		overhang_tr = data->p->max_hang < ove_tmp.mln * data->p->max_hang_ratio ? data->p->max_hang : ove_tmp.mln * data->p->max_hang_ratio;
		if (matching_bases >= data->p->m_b && ove_tmp.mln > data->p->min_ove)
		{
			// if (data->p->iter <= 1 && left_overhang < data->p->max_hang && right_overhang < data->p->max_hang)
			// if (left_overhang < overhang_tr && right_overhang < overhang_tr)
			if (overhang < overhang_tr)
			{
				ove_cl_tmp[n_tmp++] = ove_tmp;
				if (n_tmp > max_index_n) printf("[%s]2 use unlegal memory, %d-%d, something wrong with it...", __func__, n_tmp, max_index_n);
			}
			// else if (data->p->iter > 1)
			// {
			// 	ove_cl_tmp[n_tmp++] = ove_tmp;
			// }
		}
	}

	if (n_tmp == 0)
	{
		if (ove_cl_tmp != NULL) {free(ove_cl_tmp); ove_cl_tmp = NULL;}
		return 0;
	}

	qsort(ove_cl_tmp, n_tmp, sizeof(OVE_t), compare_ove_score);
	n_tmp = n_tmp > data->p->top_n ? data->p->top_n : n_tmp;

	for (i = 0; i < n_tmp; i++)
	{
		if (buf->ove_cl[q_idx].n >= buf->ove_cl[q_idx].m)
		{
			buf->ove_cl[q_idx].m = buf->ove_cl[q_idx].m << 1;
			buf->ove_cl[q_idx].ove = (OVE_t *)realloc(buf->ove_cl[q_idx].ove, buf->ove_cl[q_idx].m * sizeof(OVE_t));
			if (buf->ove_cl[q_idx].ove == NULL)
			{
				fprintf(stderr, "[%s] calloc %ldGB buf->ove_cl[q_idx].ove memory error!\n", __func__, buf->ove_cl[q_idx].m * sizeof(OVE_t) / 1024 / 1024 / 1024);
				exit(1);
			}
			printf("[%s] will not use these code, something wrong with it...",__func__);
		}
		cp_ove(&buf->ove_cl[q_idx].ove[buf->ove_cl[q_idx].n], &ove_cl_tmp[i]);
		buf->ove_cl[q_idx].ove[buf->ove_cl[q_idx].n++].tid = data->p->read_stat->iter_idx[ove_cl_tmp[i].tid];
		if (buf->ove_cl[q_idx].n > buf->ove_cl[q_idx].m) printf("[%s]buf[%d]->ove_cl[%d] memory leak..,%ld > %ld\n",__func__,pid,q_idx,buf->ove_cl[q_idx].n,buf->ove_cl[q_idx].m);
	}


	if (ove_cl_tmp != NULL) {free(ove_cl_tmp); ove_cl_tmp = NULL;}

	return n_tmp;
}

// the i is the read array's index
static void overlapping_core(void *data, int64_t i, int pid)
{
	step_query_t *step = (step_query_t *)data;
	
	finding_hits(step, i, pid);

	if (step->buf[pid]->v_n > 0)
	{
		create_graph(step->buf[pid]);
		// printf("%ld %d %d\n", i, step->buf[pid]->v_n, step->buf[pid]->max_index_n);
		get_overlap_info(step, step->buf[pid], i, pid);
	}

	return;
}

static void *mapping_pipeline(void *shared, int step, void *in)
{
	int i, j, k;
	pipeline_t *p = (pipeline_t *)shared;
	if (step == 0)
	{
		step_query_t *s = (step_query_t *)calloc(1, sizeof(step_query_t));
		s->read_info = bseq_read_map(p->fp, p->query_batch_size, &s->n_seq, &p->map_n_sta, &p->map_idx_n_sta, p->read_stat->iter_map);
		if (s->n_seq)
		{
			s->p = p;
			s->ove_alloc = s->n_seq;
			p->n_processed += s->n_seq;

			s->buf = (thread_buf_t **)calloc(p->n_threads, sizeof(thread_buf_t *));
			for (i = 0; i < p->n_threads; i++)
			{
				s->buf[i] = (thread_buf_t *)calloc(1, sizeof(thread_buf_t));
				s->buf[i]->seq_mi = (mini_t *)calloc(1, sizeof(mini_t));
				s->buf[i]->seq_mi->mi_n = 0;
				s->buf[i]->seq_mi->mi_m = p->max_len;
				s->buf[i]->seq_mi->mi_v = (MINIMIZER_t *)calloc(s->buf[i]->seq_mi->mi_m, sizeof(MINIMIZER_t));
				s->buf[i]->v_n = 0;
				s->buf[i]->v_m = p->rep_n * p->max_len;
				s->buf[i]->vertex_mr = (MR_t *)calloc(s->buf[i]->v_m, sizeof(MR_t));
				s->buf[i]->max_v_num = p->rep_n * p->max_len;
				init_graph(&s->buf[i]->Dp_graph, s->buf[i]->max_v_num);
				s->buf[i]->path = (PATH_t *)calloc(s->buf[i]->max_v_num, sizeof(PATH_t));
				s->buf[i]->max_index = (uint32_t *)calloc(s->buf[i]->max_v_num, sizeof(uint32_t));
				s->buf[i]->ove_cl = (OVE_C_t *)calloc(s->ove_alloc, sizeof(OVE_C_t));
				for (j = 0; j < s->ove_alloc; j++)
				{
					s->buf[i]->ove_cl[j].n = 0;
					s->buf[i]->ove_cl[j].m = p->top_n * 2;
					s->buf[i]->ove_cl[j].ove = (OVE_t *)calloc(s->buf[i]->ove_cl[j].m, sizeof(OVE_t));
				}
			}
			// printf("p->read_stat->iter_idx idx: %d-%d\n",p->index_idx_n_sta, s->ove_alloc + p->index_idx_n_sta - 1);
			return s;
		}
		else
		{
			if (s != NULL) {free(s); s = NULL;}
		}
	}
	else if (step == 1)
	{
		double time_beg = realtime(), cpu_beg = cputime();
		kt_for(p->n_threads, overlapping_core, in, ((step_query_t *)in)->n_seq);
		double time_end = realtime(), cpu_end = cputime();
		mapping_time += (time_end - time_beg);
		cpu_time += (cpu_end - cpu_beg);

		fprintf(stderr, "[%s] mapped %d query read...\n", __func__, ((step_query_t *)in)->n_seq);
		return in;
	}
	else if (step == 2)
	{
		uint32_t id;
		step_query_t *s = (step_query_t *)in;
		for (i = 0; i < s->n_seq; i++)
		{
			READ_t *read = &s->read_info[i];
			if (read->read_seq != NULL) {free(read->read_seq); read->read_seq = NULL;}
			if (read->read_name != NULL) {free(read->read_name); read->read_name = NULL;}
		}
		if (s->read_info != NULL) {free(s->read_info); s->read_info = NULL;}

		for (i = 0; i < p->n_threads; i++)
		{
			for (j = 0; j < s->ove_alloc; j++)
			{
				if (s->buf[i]->ove_cl[j].n == 0) continue;
				id = s->buf[i]->ove_cl[j].ove[0].qid;
				for (k = 0; k < s->buf[i]->ove_cl[j].n; k++)
				{
					if (p->ove_cl[id].n >= p->ove_cl[id].m)
					{
						p->ove_cl[id].m = p->ove_cl[id].m << 1;
						p->ove_cl[id].ove = (OVE_t *)realloc(p->ove_cl[id].ove, p->ove_cl[id].m * sizeof(OVE_t));
					}
					cp_ove(&p->ove_cl[id].ove[p->ove_cl[id].n++], &s->buf[i]->ove_cl[j].ove[k]);
					if (p->ove_cl[id].n > p->ove_cl[id].m) printf("[%s]p->ove_cl[%d] memory leak..,%ld > %ld\n",__func__,id,p->ove_cl[id].n,p->ove_cl[id].m);
				}
			}
		}

		for (i = 0; i < p->n_threads; i++)
		{
			if (s->buf[i]->seq_mi->mi_v != NULL) {free(s->buf[i]->seq_mi->mi_v);s->buf[i]->seq_mi->mi_v = NULL;}
			if (s->buf[i]->seq_mi != NULL) {free(s->buf[i]->seq_mi);s->buf[i]->seq_mi = NULL;}
			if (s->buf[i]->vertex_mr != NULL) {free(s->buf[i]->vertex_mr);s->buf[i]->vertex_mr = NULL;}
			if (s->buf[i]->path != NULL) {free(s->buf[i]->path);s->buf[i]->path = NULL;}
			if (s->buf[i]->max_index != NULL) {free(s->buf[i]->max_index);s->buf[i]->max_index = NULL;}
			free_graph(&s->buf[i]->Dp_graph);
			for (j = 0; j < s->ove_alloc; j++)
			{
				if (s->buf[i]->ove_cl[j].ove != NULL) {free(s->buf[i]->ove_cl[j].ove);s->buf[i]->ove_cl[j].ove = NULL;}
			}
			if (s->buf[i]->ove_cl != NULL) {free(s->buf[i]->ove_cl);s->buf[i]->ove_cl = NULL;}
			if (s->buf[i] != NULL) {free(s->buf[i]);s->buf[i] = NULL;}
		}
		if (s->buf != NULL) {free(s->buf); s->buf = NULL;}
		if (s != NULL) {free(s); s = NULL;}
	}
	return 0;
}

uint32_t estimate_read_portion_coverage(pipeline_t *p, read_cov_t *read_cov, uint32_t split_len, double *all_ave_cov)
{
	uint32_t seed_rid;
	int32_t cov_q_sta, cov_q_end, cov_q_sta_pos, cov_q_end_pos;
	int32_t cov_sta, cov_end; // cutting region out of seed read
	cov_data_t last_read, cur_read;
	uint32_t *visited_rid = (uint32_t *)calloc(p->n_total_read, sizeof(uint32_t));
	uint32_t size_tr = split_len / 2, last_size, sum_cov, cov_nn = p->max_len / split_len + 1;
	double ave_cov;
	uint32_t *cov_n = (uint32_t *)calloc(cov_nn, sizeof(uint32_t));
    rStack stk;
    kv_init(stk);
	uint32_t cov_limit = 256, cov_acc = 0;

	// step 1: estimate coverage of seed read
	for(uint32_t i = 0; i < p->n_total_read; ++i)
	{
		if (p->read_stat->is_idx[i] != 1) continue;
    	kv_init(stk);
		seed_rid = i;
		printf("seed read %d(%d)\t", seed_rid, read_cov[seed_rid].rid);
		cov_q_sta_pos = cov_q_end_pos = 0;
		last_read.rid = seed_rid;
		last_read.sta_pos = last_read.end_pos = 0;
		last_read.ts = last_read.te = 0;
		last_read.rev = 0;
		kv_push(cov_data_t, stk, last_read);
		visited_rid[last_read.rid] = 1;
		cov_acc = 0;
		while (kv_size(stk) > 0)
		{
			last_read = kv_pop(stk);
			// if (seed_rid == 27) printf("pop read %d\n", last_read.rid);
			for(uint32_t j = 0; j < p->ove_cl[last_read.rid].nn; ++ j)
			{
				if (visited_rid[p->ove_cl[last_read.rid].ove[j].tid] == 1)	continue;
				// if (seed_rid == 27) printf("%u\t%u\t%u\t%u\t%c\t%u\t%u\t%u\t%u\t%d\t%d\t%.2f\n", p->ove_cl[last_read.rid].ove[j].qid, p->ove_cl[last_read.rid].ove[j].ql, p->ove_cl[last_read.rid].ove[j].qs, p->ove_cl[last_read.rid].ove[j].qe, "+-"[p->ove_cl[last_read.rid].ove[j].rev], p->ove_cl[last_read.rid].ove[j].tid, p->ove_cl[last_read.rid].ove[j].tl, p->ove_cl[last_read.rid].ove[j].ts, p->ove_cl[last_read.rid].ove[j].te, p->ove_cl[last_read.rid].ove[j].mbp, p->ove_cl[last_read.rid].ove[j].mln, p->ove_cl[last_read.rid].ove[j].score);
				// if (seed_rid == 27) printf("rev: %d\n", last_read.rev);
				if (last_read.rev == 1)
				{
					cov_q_sta_pos = last_read.sta_pos - (p->ove_cl[last_read.rid].ove[j].qe - last_read.te);
					cov_q_end_pos = last_read.end_pos - (p->ove_cl[last_read.rid].ove[j].qs - last_read.ts);
					// if (seed_rid == 27) printf("sta = %d - (%d - %d) = %d\n", last_read.sta_pos, p->ove_cl[last_read.rid].ove[j].qe, last_read.te, cov_q_sta_pos);
					// if (seed_rid == 27) printf("end = %d - (%d - %d) = %d\n", last_read.end_pos, p->ove_cl[last_read.rid].ove[j].qs, last_read.ts, cov_q_end_pos);
				}
				else
				{
					cov_q_sta_pos = last_read.sta_pos + p->ove_cl[last_read.rid].ove[j].qs - last_read.ts;
					cov_q_end_pos = last_read.end_pos + p->ove_cl[last_read.rid].ove[j].qe - last_read.te;
					// if (seed_rid == 27) printf("sta = %d + %d - %d = %d\n", last_read.sta_pos, p->ove_cl[last_read.rid].ove[j].qs, last_read.ts, cov_q_sta_pos);
					// if (seed_rid == 27) printf("end = %d + %d - %d = %d\n", last_read.end_pos, p->ove_cl[last_read.rid].ove[j].qe, last_read.te, cov_q_end_pos);
				}

				// if (seed_rid == 27) printf("q[%d-%d]->q[%d-%d]\n", last_read.sta_pos, last_read.end_pos, cov_q_sta_pos, cov_q_end_pos);
				if (cov_q_sta_pos > (int32_t)p->read_stat->rlen[seed_rid] || cov_q_end_pos < 0)	continue;
				cov_sta = cov_q_sta_pos < 0 ? 0 : cov_q_sta_pos;
				cov_end = cov_q_end_pos > (int32_t)p->read_stat->rlen[seed_rid] ? (int32_t)p->read_stat->rlen[seed_rid] : cov_q_end_pos;
				last_size = (p->read_stat->rlen[seed_rid] % split_len);
				last_size = last_size > (split_len / 2) ? last_size : 0;
				cov_q_sta = cov_sta % split_len < size_tr ? cov_sta / split_len : cov_sta / split_len + 1;
				if (cov_end / split_len == read_cov[seed_rid].n - 1) // the last region of seed read < split_len bp
					cov_q_end = cov_end % split_len > 0 ? read_cov[seed_rid].n - 1 : read_cov[seed_rid].n - 2;
				else
					cov_q_end = cov_end % split_len > size_tr ? cov_end / split_len : cov_end / split_len - 1;

				// if (seed_rid == 27) printf("q[%d-%d]\n", cov_q_sta, cov_q_end);
				// if (cov_q_sta > cov_q_end) continue;
				for(int32_t k = cov_q_sta; k <= cov_q_end; ++k)
				{
					read_cov[seed_rid].cov[k] += 1;
					cov_acc = cov_acc > read_cov[seed_rid].cov[k] ? cov_acc : read_cov[seed_rid].cov[k];
				}
				// if (seed_rid == 27) for(uint32_t k = 0; k < read_cov[seed_rid].n; ++k) printf("%d\t",read_cov[seed_rid].cov[k]);
				// if (seed_rid == 27) printf("\n");

				if (p->read_stat->is_idx[p->ove_cl[last_read.rid].ove[j].tid] == 1)
				{
					cur_read.rid = p->ove_cl[last_read.rid].ove[j].tid;
					cur_read.sta_pos = cov_q_sta_pos;
					cur_read.end_pos = cov_q_end_pos;
					cur_read.ts = p->ove_cl[last_read.rid].ove[j].ts > p->ove_cl[last_read.rid].ove[j].te ? p->ove_cl[last_read.rid].ove[j].te : p->ove_cl[last_read.rid].ove[j].ts;
					cur_read.te = p->ove_cl[last_read.rid].ove[j].ts > p->ove_cl[last_read.rid].ove[j].te ? p->ove_cl[last_read.rid].ove[j].ts : p->ove_cl[last_read.rid].ove[j].te;
					cur_read.rev = last_read.rev ^ p->ove_cl[last_read.rid].ove[j].rev;
					kv_push(cov_data_t, stk, cur_read);
				}
				visited_rid[p->ove_cl[last_read.rid].ove[j].tid] = 1;
			}
			if (cov_acc > cov_limit) break;
		}
		sum_cov = 0;
		for(uint32_t k = 0; k < read_cov[seed_rid].n; ++k)
		{
			sum_cov += read_cov[seed_rid].cov[k];
			printf("%d\t",read_cov[seed_rid].cov[k]);
		}
		read_cov[seed_rid].ave_cov = sum_cov / (double)read_cov[seed_rid].n;
		printf("ave %.3f\n", read_cov[seed_rid].ave_cov);
		*all_ave_cov += read_cov[seed_rid].ave_cov;
		memset(visited_rid, 0, p->n_total_read * sizeof(uint32_t));
	}

	int32_t cov_t_sta, cov_t_end, cov_t_sta_pos, cov_t_end_pos;
	// step 2: estimate coverage of other read
	for(uint32_t i = 0; i < p->n_total_read; ++i)
	{
		if (p->read_stat->is_idx[i] == 1) continue;
		for (uint32_t j = 0; j < cov_nn; ++j) cov_n[j] = 0;
		for (uint32_t j = 0; j < p->ove_cl[i].nn; ++j)
		{
			// if (i == 10045) printf("%u\t%u\t%u\t%u\t%c\t%u\t%u\t%u\t%u\t%d\t%d\t%.2f\n", p->ove_cl[i].ove[j].qid, p->ove_cl[i].ove[j].ql, p->ove_cl[i].ove[j].qs, p->ove_cl[i].ove[j].qe, "+-"[p->ove_cl[i].ove[j].rev], p->ove_cl[i].ove[j].tid, p->ove_cl[i].ove[j].tl, p->ove_cl[i].ove[j].ts, p->ove_cl[i].ove[j].te, p->ove_cl[i].ove[j].mbp, p->ove_cl[i].ove[j].mln, p->ove_cl[i].ove[j].score);
			seed_rid = p->ove_cl[i].ove[j].tid;
			cov_t_sta_pos = p->ove_cl[i].ove[j].rev == 0 ? p->ove_cl[i].ove[j].ts : p->ove_cl[i].ove[j].te;
			cov_t_end_pos = p->ove_cl[i].ove[j].rev == 0 ? p->ove_cl[i].ove[j].te : p->ove_cl[i].ove[j].ts;
			cov_t_sta = cov_t_sta_pos / split_len;
			cov_t_end = cov_t_end_pos / split_len;
			// if (i == 10045) printf("t[%d-%d]->t[%d-%d]\n", cov_t_sta_pos, cov_t_end_pos, cov_t_sta, cov_t_end);
			sum_cov = 0;
			for (int32_t k = cov_t_sta; k <= cov_t_end; ++k)
			{
				sum_cov += read_cov[seed_rid].cov[k];
				// if (i == 10045) printf("%d\t", read_cov[seed_rid].cov[k]);
			}
			ave_cov = sum_cov / (double)(cov_t_end - cov_t_sta + 1);
			// if (i == 10045) printf("ave %.3f\n", ave_cov);

			last_size = (p->read_stat->rlen[seed_rid] % split_len);
			cov_q_sta_pos = p->ove_cl[i].ove[j].rev == 0 ? p->ove_cl[i].ove[j].qs - cov_t_sta_pos : p->ove_cl[i].ove[j].qs - (p->read_stat->rlen[seed_rid] - cov_t_end_pos);
			cov_q_end_pos = p->ove_cl[i].ove[j].rev == 0 ? p->ove_cl[i].ove[j].qe + (p->read_stat->rlen[seed_rid] - cov_t_end_pos) : p->ove_cl[i].ove[j].qe + cov_t_sta_pos;
			cov_sta = cov_q_sta_pos < 0 ? 0 : cov_q_sta_pos;
			cov_end = cov_q_end_pos > (int32_t)p->read_stat->rlen[i] ? (int32_t)p->read_stat->rlen[i] : cov_q_end_pos;
			cov_q_sta = cov_sta / split_len;
			cov_q_end = cov_end / split_len;
			// if (i == 10045) printf("q[%d-%d]->q[%d-%d]->q[%d-%d]\n", cov_q_sta_pos, cov_q_end_pos, cov_sta, cov_end, cov_q_sta, cov_q_end);
			for (int32_t k = cov_q_sta; k <= cov_q_end; ++k)
			{
				read_cov[i].cov[k] += ave_cov;
				cov_n[k] += 1;
				// if (i == 10045) printf("%d\t", read_cov[i].cov[k]);
			}
			// if (i == 10045) printf("ave %.3f\n", ave_cov);
		}
		sum_cov = 0;
		printf("other read %d(%d)\t", i, read_cov[i].rid);
		for (uint32_t k = 0; k < read_cov[i].n; ++k)
		{
			if (cov_n[k] == 0)
			{
				printf("%d\t", read_cov[i].cov[k]);
				continue;
			}
			read_cov[i].cov[k] = read_cov[i].cov[k] / cov_n[k];
			sum_cov += read_cov[i].cov[k];
			printf("%d\t", read_cov[i].cov[k]);
		}
		read_cov[i].ave_cov = sum_cov / (double)read_cov[i].n;
		printf("ave %.3f\n", read_cov[i].ave_cov);
		*all_ave_cov += read_cov[i].ave_cov;
	}

	kv_destroy(stk);
	if (visited_rid != NULL) {free(visited_rid); visited_rid = NULL;}
	if (cov_n != NULL) {free(cov_n); cov_n = NULL;}

	return 0;
}

int is_redundant_ove(OVE_C_t *seed_ove_cl, uint32_t r_i)
{
	for(uint32_t i = 0; i < seed_ove_cl->n; ++i)
	{
		if (seed_ove_cl->ove[i].tid == r_i)	return 1;
	}
	return 0;
}

int mark_read_not_idx(uint32_t rid, OVE_C_t *ove_cl, read_stat_t *read_stat)
{
	uint32_t tid;
	for(uint32_t i = 0; i < ove_cl[rid].n; ++i)
	{
		tid = ove_cl[rid].ove[i].tid;
		read_stat->not_idx[tid] = 1;
		printf("not idx read %d = 1\n",tid);
	}
	return 0;
}

int choose_idx_read_for_next_iteration(pipeline_t *pl, read_cov_t *read_cov, double all_ave_cov)
{
	int ret, tmp;
	uint32_t buf_size = part_size;
	uint32_t *buf = (uint32_t *)calloc(buf_size, sizeof(uint32_t));
	uint32_t buf_pos, buf_cov, ele;
	double buf_ave_cov;

	pl->read_stat->iter_idx_n = 0;
	pl->read_stat->iter_map_n = 0;
	memset(pl->read_stat->not_idx, 0, pl->n_total_read * sizeof(uint32_t));

	// identify region(>5000bp) with not enough coverage
	for (int r_i = 0; r_i < pl->n_total_read; r_i++)
	{
		srand(r_i * (pl->iter + 1));
		// printf("new read %d, %d block\n", r_i, read_cov[r_i].n);
		ret = 0;
		buf_pos = 0;
		buf_cov = 0;
		for (int k = 0; k < read_cov[r_i].n; k++)
		{
			ele = buf[buf_pos];
			buf[buf_pos] = read_cov[r_i].cov[k];
			buf_cov += buf[buf_pos++];
			if (k >= buf_size)
			{
				buf_cov -= ele;
				buf_ave_cov = (double)buf_cov / buf_size;
				// printf("%d %f\n", k, buf_ave_cov);
				if (buf_ave_cov <= all_ave_cov * cov_ratio)	{ret = 1; break;}
			}
			if (buf_pos == buf_size) buf_pos = 0;
		}
		if (read_cov[r_i].n <= buf_size)
		{
			buf_ave_cov = buf_cov / read_cov[r_i].n;
			// printf("%d %f\n", read_cov[r_i].n, buf_ave_cov);
			if (buf_ave_cov <= all_ave_cov * cov_ratio)	ret = 1;
		}
		if (ret == 1)
		{
			pl->read_stat->iter_map[pl->read_stat->iter_map_n++] = read_cov[r_i].rid;
			tmp = rand() % 100;
			if (tmp < X_read && pl->read_stat->not_idx[r_i] != 1)
			{
				pl->read_stat->iter_idx[pl->read_stat->iter_idx_n++] = read_cov[r_i].rid;
				// mark_read_not_idx(read_cov[r_i].rid, pl->ove_cl, pl->read_stat);
				pl->read_stat->is_idx[r_i] = 1;
			}
		}
		// if (r_i == 380) break;
	}

	for (int r_i = 0; r_i < pl->read_stat->iter_idx_n; r_i++)
	{
		printf("idx read %d\t\n",pl->read_stat->iter_idx[r_i]);
		for(uint32_t k = 0; k < read_cov[pl->read_stat->iter_idx[r_i]].n; ++k)
		{
			printf("%d\t",read_cov[pl->read_stat->iter_idx[r_i]].cov[k]);
		}
		printf("ave %.3f\n", read_cov[pl->read_stat->iter_idx[r_i]].ave_cov);
	}

	for (int r_i = 0; r_i < pl->read_stat->iter_map_n; r_i++)
	{
		printf("map read %d\t\n",pl->read_stat->iter_map[r_i]);
		for(uint32_t k = 0; k < read_cov[pl->read_stat->iter_map[r_i]].n; ++k)
		{
			printf("%d\t",read_cov[pl->read_stat->iter_map[r_i]].cov[k]);
		}
		printf("ave %.3f\n", read_cov[pl->read_stat->iter_map[r_i]].ave_cov);
	}
	fprintf(stderr, "[To next iteration] %d reads do not have enough coverage( < %f), index %d read...\n", pl->read_stat->iter_map_n, all_ave_cov * cov_ratio, pl->read_stat->iter_idx_n);

	if (buf != NULL) {free(buf); buf = NULL;}
	return 0;
}

uint32_t finding_overlapping(const char *read_fastq, const char *index_fastq, const char *temp_ove_dir, const char *temp_ove_link_dir)
{
	uint32_t split_len = 1000;
	double all_ave_cov = 0;
	read_cov_t *read_cov;

	// step 1: sorting reads' length
	bseq_file_t *bf = NULL;
	bf = bseq_open(read_fastq);
	if (bf == NULL)
	{
		fprintf(stderr, "[Warning] Wrong input file route or name:%s \n", read_fastq);
		exit(1);
	}

	pipeline_t pl;
	pl.x_len = 0;
	pl.max_len = 0;
	pl.batch_i = 0;
	pl.n_processed = 0;
	pl.n_total_read = 0;
	pl.index_batch_size = batch_size_read;
	pl.fp = bf;
	pl.read_stat = (read_stat_t *)calloc(1, sizeof(read_stat_t));
	pl.read_stat->rlen = NULL;
	pl.read_stat->is_idx = NULL;
	pl.read_stat->not_idx = NULL;
	pl.read_stat->iter_idx = NULL;
	pl.read_stat->iter_idx_n = 0;
	pl.read_stat->iter_map = NULL;
	pl.read_stat->iter_map_n = 0;

	kt_pipeline(thread_n < 2 ? thread_n : 2, sort_read_length_pipeline, &pl, 2);

	pl.x_len = pl.x_len / pl.batch_i;
	pl.n_total_read = pl.n_processed;
	if (bf != NULL) bseq_close(bf);

	fprintf(stderr, "[Index] Sorted %d read length, indexing read > %d..., max read length %d\n", pl.n_total_read, pl.x_len, pl.max_len);
	fprintf(stderr, "[%s] Real time:%.6f sec\tCPU:%.6f sec\tMemory peak:%.6f GB\n", __func__, realtime() - realtime0, cputime(), peak_memory() / 1024.0 / 1024.0);

	// step 2: iteratively indexing and mapping
	pl.iter = 0;
	pl.n_threads = thread_n;
	pl.query_batch_size = 1000000000;
	pl.m_b = m_b;
	pl.min_ove = min_ove;
	pl.max_hang = max_hang;
	pl.max_hang_ratio = max_hang_ratio;
	pl.top_n = top_n;
	pl.ove_cl = (OVE_C_t *)calloc(pl.n_total_read, sizeof(OVE_C_t));
	if (pl.ove_cl == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB pl.ove_cl memory error!\n", __func__, pl.n_total_read * sizeof(OVE_C_t) / 1024 / 1024 / 1024);
		exit(1);
	}
	for (uint32_t i = 0; i < pl.n_total_read; i++)
	{
		pl.ove_cl[i].n = 0;
		pl.ove_cl[i].m = pl.top_n;
		pl.ove_cl[i].ove = (OVE_t *)calloc(pl.ove_cl[i].m, sizeof(OVE_t));
	}
	fprintf(stderr, "[%s] calloc %fGB pl.ove_cl[i].ove memory...\n", __func__, pl.n_total_read * (pl.top_n * sizeof(OVE_t) + 5 * sizeof(uint32_t)) / 1024.0 / 1024.0 / 1024.0);

	pl.read_stat->is_idx = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	pl.read_stat->not_idx = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	pl.read_stat->iter_idx = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	pl.read_stat->iter_map = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	for (int r_i = 0; r_i < pl.n_total_read; ++r_i)
	{
		if (pl.read_stat->rlen[r_i] >= pl.x_len)
		{
			pl.read_stat->is_idx[r_i] = 1;
			pl.read_stat->iter_idx[pl.read_stat->iter_idx_n++] = r_i; // mark seed read for the first iteration
		}
		pl.read_stat->iter_map[pl.read_stat->iter_map_n++] = r_i; // mark seed read for the first iteration
	}
	// for (int r_i = 0; r_i < pl.read_stat->iter_idx_n; r_i++) printf("idx read\t%d\t\n",pl.read_stat->iter_idx[r_i]);
	// for (int r_i = 0; r_i < pl.read_stat->iter_map_n; r_i++) printf("map read\t%d\t\n",pl.read_stat->iter_map[r_i]);

	read_cov = (read_cov_t *)calloc(pl.n_total_read, sizeof(read_cov_t));
	for (int r_i = 0; r_i < pl.n_total_read; ++r_i)
	{
		read_cov[r_i].rid = r_i;
		read_cov[r_i].ave_cov = 0.;
		read_cov[r_i].n = pl.read_stat->rlen[r_i] / split_len + 1;
		read_cov[r_i].cov = (uint32_t *)calloc(read_cov[r_i].n, sizeof(uint32_t));
	}

	bseq_file_t *fp = NULL;

	int ex_iter = 2;
	while(pl.iter < ex_iter)
	{
		pl.index_n_sta = 0;
		pl.index_idx_n_sta = 0;
		pl.index_idx_n_end = 0;
		fp = bseq_open(index_fastq);

		for(;;)
		{
			// step 2.1: indexing reads in read_stat->iter_idx
			mini_idx_t *mi = NULL;
			pl.index_idx_n_sta = pl.index_idx_n_end;
			mi = indexing_input_read(fp, &pl.rep_n, &pl.index_n_sta, &pl.index_idx_n_end, pl.read_stat);
			// fprintf(stderr, "index idx sta %d, end %d\n",pl.index_idx_n_sta,pl.index_idx_n_end);
			if (mi->mm.mi_n == 0)
			{
				if (mi->mm.mi_v != NULL) {free(mi->mm.mi_v); mi->mm.mi_v = NULL;}
				if (mi->mi_count != NULL) {free(mi->mi_count); mi->mi_count = NULL;}
				if (mi != NULL) {free(mi); mi = NULL;}
				break;
			}

			// step 2.2: mapping reads to indexed reads
			bf = bseq_open(read_fastq);
			pl.mi = mi;
			pl.fp = bf;
			pl.n_processed = 0;
			pl.map_n_sta = 0;
			pl.map_idx_n_sta = 0;
			pl.map_idx_n_end = 0;

			kt_pipeline(thread_n < 3 ? thread_n : 3, mapping_pipeline, &pl, 3);
			fprintf(stderr, "[Overlap] Mapped %d query read...\n", pl.n_processed);
			fprintf(stderr, "[%s] Real time:%.6f sec\tCPU:%.6f sec\tMemory peak:%.6f GB\n", __func__, realtime() - realtime0, cputime(), peak_memory() / 1024.0 / 1024.0);

			if (bf != NULL)	bseq_close(bf);

			if (mi->mm.mi_v != NULL) {free(mi->mm.mi_v); mi->mm.mi_v = NULL;}
			if (mi->mi_count != NULL) {free(mi->mi_count); mi->mi_count = NULL;}
			if (mi != NULL) {free(mi); mi = NULL;}
		}

		if (fp != NULL) bseq_close(fp);
		
		// symm the ove_cl
		for (int r_i = 0; r_i < pl.n_total_read; ++r_i) pl.ove_cl[r_i].nn = pl.ove_cl[r_i].n;
		
		for (int r_i = 0; r_i < pl.n_total_read; ++r_i)
		{
			if (pl.ove_cl[r_i].n == 0) continue;
			for (int o_i = 0; o_i < pl.ove_cl[r_i].n; ++o_i)
			{
				int seed_rid = pl.ove_cl[r_i].ove[o_i].tid;
				if (pl.ove_cl[seed_rid].nn >= pl.ove_cl[seed_rid].m)
				{
					pl.ove_cl[seed_rid].m = pl.ove_cl[seed_rid].m << 1;
					pl.ove_cl[seed_rid].ove = (OVE_t *)realloc(pl.ove_cl[seed_rid].ove, pl.ove_cl[seed_rid].m * sizeof(OVE_t));
				}

				if (pl.read_stat->is_idx[r_i] == 1 && is_redundant_ove(&pl.ove_cl[seed_rid], r_i) == 1) continue;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn].qid = pl.ove_cl[r_i].ove[o_i].tid;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn].ql = pl.ove_cl[r_i].ove[o_i].tl;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn].qs = pl.ove_cl[r_i].ove[o_i].ts > pl.ove_cl[r_i].ove[o_i].te ? pl.ove_cl[r_i].ove[o_i].te : pl.ove_cl[r_i].ove[o_i].ts;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn].qe = pl.ove_cl[r_i].ove[o_i].ts > pl.ove_cl[r_i].ove[o_i].te ? pl.ove_cl[r_i].ove[o_i].ts : pl.ove_cl[r_i].ove[o_i].te;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn].tid = pl.ove_cl[r_i].ove[o_i].qid;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn].tl = pl.ove_cl[r_i].ove[o_i].ql;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn].ts = pl.ove_cl[r_i].ove[o_i].ts > pl.ove_cl[r_i].ove[o_i].te ? pl.ove_cl[r_i].ove[o_i].qe : pl.ove_cl[r_i].ove[o_i].qs;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn].te = pl.ove_cl[r_i].ove[o_i].ts > pl.ove_cl[r_i].ove[o_i].te ? pl.ove_cl[r_i].ove[o_i].qs : pl.ove_cl[r_i].ove[o_i].qe;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn].rev = pl.ove_cl[r_i].ove[o_i].rev;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn].mbp = pl.ove_cl[r_i].ove[o_i].mbp;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn].mln = pl.ove_cl[r_i].ove[o_i].mln;
				pl.ove_cl[seed_rid].ove[pl.ove_cl[seed_rid].nn++].score = pl.ove_cl[r_i].ove[o_i].score;
			}
		}

		// generate read portion coverage statistics from overlap collection
		// re initial array read_cov[]
		for (int r_i = 0; r_i < pl.n_total_read; ++r_i) {read_cov[r_i].ave_cov = 0.; memset(read_cov[r_i].cov, 0, read_cov[r_i].n * sizeof(uint32_t));}
		estimate_read_portion_coverage(&pl, read_cov, split_len, &all_ave_cov);
		all_ave_cov = all_ave_cov / pl.n_total_read;
		fprintf(stderr, "[Current iteration] average coverage of all reads: %f...\n",all_ave_cov);
		// if (pl.iter < ex_iter - 1)
		// {
			choose_idx_read_for_next_iteration(&pl, read_cov, all_ave_cov);
		// }

		if (pl.read_stat->iter_idx_n == 0) break;
		
		pl.iter++;
	}

	fprintf(stderr, "[Statistic] Mapping time:%.6f sec\tCPU:%.6f sec\tMemory peak:%.6f GB\n", mapping_time, cpu_time, peak_memory() / 1024.0 / 1024.0);

	// step 4: output
	FILE *fp_ove = NULL;
	FILE *fp_ove_link = NULL;

	fp_ove_link = fopen(temp_ove_link_dir, "w");
	if (fp_ove_link == NULL)
	{
		fprintf(stderr, "[%s Wrong] Failed to open file %s!!!\n", __func__, temp_ove_link_dir);
		exit(1);
	}

	fp_ove = fopen(temp_ove_dir, "w");
	if (fp_ove == NULL)
	{
		fprintf(stderr, "[%s Wrong] Failed to open file %s!!!\n",  __func__, temp_ove_dir);
		exit(1);
	}

	uint32_t total_index_num = 0;
	for (int r_i = 0; r_i < pl.n_total_read; ++r_i)
	{
		if (pl.read_stat->is_idx[r_i] == 1) total_index_num += 1;
		if (pl.ove_cl[r_i].n == 0) continue;
		for (int o_i = 0; o_i < pl.ove_cl[r_i].n; ++o_i)
		{
			// if (pl.ove_cl[r_i].ove[o_i].ql < pl.x_len)
			// {
			// 	fprintf(fp_ove_link, "%u\t%u\t%u\t%u\t%c\t", pl.ove_cl[r_i].ove[o_i].qid, pl.ove_cl[r_i].ove[o_i].ql, pl.ove_cl[r_i].ove[o_i].qs, pl.ove_cl[r_i].ove[o_i].qe, "+-"[pl.ove_cl[r_i].ove[o_i].rev]);
			// 	fprintf(fp_ove_link, "%u\t%u\t%u\t%u\t", pl.ove_cl[r_i].ove[o_i].tid, pl.ove_cl[r_i].ove[o_i].tl, pl.ove_cl[r_i].ove[o_i].ts, pl.ove_cl[r_i].ove[o_i].te);
			// 	fprintf(fp_ove_link, "%d\t%d\t%.2f\n", pl.ove_cl[r_i].ove[o_i].mbp, pl.ove_cl[r_i].ove[o_i].mln, pl.ove_cl[r_i].ove[o_i].score);
			// }
			// if (pl.ove_cl[r_i].ove[o_i].ql >= pl.x_len)
			// {
				fprintf(fp_ove, "%u\t%u\t%u\t%u\t%c\t", pl.ove_cl[r_i].ove[o_i].qid, pl.ove_cl[r_i].ove[o_i].ql, pl.ove_cl[r_i].ove[o_i].qs, pl.ove_cl[r_i].ove[o_i].qe, "+-"[pl.ove_cl[r_i].ove[o_i].rev]);
				fprintf(fp_ove, "%u\t%u\t%u\t%u\t", pl.ove_cl[r_i].ove[o_i].tid, pl.ove_cl[r_i].ove[o_i].tl, pl.ove_cl[r_i].ove[o_i].ts, pl.ove_cl[r_i].ove[o_i].te);
				fprintf(fp_ove, "%d\t%d\t%.2f\ti:%c\n", pl.ove_cl[r_i].ove[o_i].mbp, pl.ove_cl[r_i].ove[o_i].mln, pl.ove_cl[r_i].ove[o_i].score, "mi"[pl.read_stat->is_idx[r_i]]);
			// }
		}
		for (int o_i = 0; o_i < pl.ove_cl[r_i].nn; ++o_i)
		{
				fprintf(fp_ove_link, "%u\t%u\t%u\t%u\t%c\t", pl.ove_cl[r_i].ove[o_i].qid, pl.ove_cl[r_i].ove[o_i].ql, pl.ove_cl[r_i].ove[o_i].qs, pl.ove_cl[r_i].ove[o_i].qe, "+-"[pl.ove_cl[r_i].ove[o_i].rev]);
				fprintf(fp_ove_link, "%u\t%u\t%u\t%u\t", pl.ove_cl[r_i].ove[o_i].tid, pl.ove_cl[r_i].ove[o_i].tl, pl.ove_cl[r_i].ove[o_i].ts, pl.ove_cl[r_i].ove[o_i].te);
				fprintf(fp_ove_link, "%d\t%d\t%.2f\ti:%c\n", pl.ove_cl[r_i].ove[o_i].mbp, pl.ove_cl[r_i].ove[o_i].mln, pl.ove_cl[r_i].ove[o_i].score, "mi"[pl.read_stat->is_idx[r_i]]);
		}
	}

	fclose(fp_ove_link);
	fclose(fp_ove);

	fprintf(stderr, "[Result] Total indexed read num: %d...\n", total_index_num);

	for (int i = 0; i < pl.n_total_read; i++)
	{
		if (pl.ove_cl[i].ove != NULL) {free(pl.ove_cl[i].ove); pl.ove_cl[i].ove = NULL;}
	}
	if (pl.ove_cl != NULL) {free(pl.ove_cl); pl.ove_cl = NULL;}
	
	for (int r_i = 0; r_i < pl.n_total_read; r_i++)
	{
		if (read_cov[r_i].cov != NULL) {free(read_cov[r_i].cov); read_cov[r_i].cov = NULL;}
	}
	if (read_cov != NULL) {free(read_cov); read_cov = NULL;}

	if (pl.read_stat->rlen != NULL) {free(pl.read_stat->rlen); pl.read_stat->rlen = NULL;}
	if (pl.read_stat->is_idx != NULL) {free(pl.read_stat->is_idx); pl.read_stat->is_idx = NULL;}
	if (pl.read_stat->not_idx != NULL) {free(pl.read_stat->not_idx); pl.read_stat->not_idx = NULL;}
	if (pl.read_stat->iter_idx != NULL) {free(pl.read_stat->iter_idx); pl.read_stat->iter_idx = NULL;}
	if (pl.read_stat->iter_map != NULL) {free(pl.read_stat->iter_map); pl.read_stat->iter_map = NULL;}
	if (pl.read_stat != NULL) {free(pl.read_stat); pl.read_stat = NULL;}

	return 0;
}