#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
// #include <string.h>
#include "bseq.h"
#include "kseq.h"
#include "overlapping.h"
#include "bit_operation.h"

#define FILTER_LEN 1000

KSEQ_INIT(gzFile, gzread)

struct bseq_file_s {
	gzFile fp;
	kseq_t *ks;
};

bseq_file_t *bseq_open(const char *fn)
{
	bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
	if (f == 0) return 0;
	fp = (bseq_file_t*)calloc(1, sizeof(bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void bseq_close(bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

// int trans_read_bit(READ_t *read_info, char *read_seq, uint32_t read_length, uint64_t *total_mem)
// {
//     uint32_t r_i = 0;
//     uint8_t c_tmp = 0;
//     char tmp_char;

//     read_info->read_bit = (uint64_t *)calloc((((read_length - 1) >> 5) + 1) << 3, sizeof(uint64_t));
//     if (read_info->read_bit == NULL)
//     {
//         fprintf(fp_debug, "calloc read_bit[0] memory error\n");
//         exit(1);
//     }

//     *total_mem += sizeof(uint64_t) * (((read_length - 1) >> 5) + 1) << 3;

//     r_i = 0;
//     while (read_seq[r_i])
//     {
//         tmp_char = read_seq[r_i];
//         if (tmp_char == 'N')
//         {
//             tmp_char = "ACGT"[rand() % 4]; //random
//         }
//         c_tmp = seq_nt4_table[(uint8_t)tmp_char];

//         read_info->read_bit[r_i >> 5] |= (((uint64_t)c_tmp) << ((31 - (r_i & 0X1f)) << 1));

//         r_i++;
//     }
//     return 0;
// }

uint32_t *bseq_read_l(bseq_file_t *fp, uint32_t batch, uint32_t *_n)
{
    uint32_t *r_len = NULL;
    uint32_t seqii = 0, m = 0;
    int64_t kr1 = 1;
    kseq_t *seq = fp->ks;

    for (seqii = 0; (seqii < batch) && ((kr1 = kseq_read(seq)) > 0); seqii++)
    {
        if (seqii >= m)
        {
            m = m == 0 ? 256 : m << 1;
            r_len = (uint32_t *)realloc(r_len, m * sizeof(uint32_t));
            if (r_len == NULL)
            {
                fprintf(stderr, "[%s] calloc %ldGB r_len memory error!\n", __func__, m * sizeof(uint32_t) / 1024 /1024 / 1024);
                exit(1);
            }
        }
        r_len[seqii] = seq->seq.l;
    }

    *_n = seqii;
    return r_len;
}

READ_t *bseq_read(bseq_file_t *fp, uint32_t batch, uint32_t *_n)
{
    READ_t *read_info = NULL;
    uint32_t idx = 0, base = 0, m = 0;
    int64_t kr1 = 1;
    kseq_t *seq1 = fp->ks;

    for (; (base < batch) && ((kr1 = kseq_read(seq1)) > 0); )
    {
        if (idx >= m)
        {
            m = m == 0 ? 256 : m << 1;
            read_info = (READ_t *)realloc(read_info, m * sizeof(READ_t));
            if (read_info == NULL)
            {
                fprintf(stderr, "[%s] calloc %ldGB read_info memory error!\n", __func__, m * sizeof(READ_t) / 1024 /1024 / 1024);
                exit(1);
            }
        }
        read_info[idx].read_name = strdup(seq1->name.s);
        read_info[idx].read_length = seq1->seq.l;
        read_info[idx].read_seq = strdup(seq1->seq.s);
        base += read_info[idx++].read_length;
    }

    *_n = idx;
    return read_info;
}

// the rid(seqii) is the same of all kind of reads
READ_t *bseq_read_idx(bseq_file_t *fp, uint32_t batch, uint32_t *_n, uint32_t *_sta, uint32_t *_sta_idx, uint32_t *iter_idx)
{
    READ_t *read_info = NULL;
    uint32_t idx = 0, base = 0, m = 0, seqi = *_sta, seqidx = *_sta_idx;
    int64_t kr1 = 1;
    kseq_t *seq1 = fp->ks;
    // fprintf(stderr, "0 seqi %d, seqidx %d, map to read id %d\n",seqi,seqidx,iter_idx[seqidx]);
    for (; (base < batch) && ((kr1 = kseq_read(seq1)) > 0); seqi++)
    {
        if (seqi == iter_idx[seqidx])
        {
            if (idx >= m)
            {
                m = m == 0 ? 256 : m << 1;
                read_info = (READ_t *)realloc(read_info, m * sizeof(READ_t));
                if (read_info == NULL)
                {
                    fprintf(stderr, "[%s] calloc %ldGB read_info memory error!\n", __func__, m * sizeof(READ_t) / 1024 /1024 / 1024);
                    exit(1);
                }
            }
            read_info[idx].rid = seqi;
            read_info[idx].read_name = strdup(seq1->name.s);
            read_info[idx].read_length = seq1->seq.l;
            read_info[idx].read_seq = strdup(seq1->seq.s);
            // is_idx[seqi] = 1;
            base += read_info[idx++].read_length;
            seqidx++;
        }
    }
    *_n = idx;
    *_sta = seqi;
    *_sta_idx = seqidx;
    // fprintf(stderr, "1 seqi %d, seqidx %d, map to read id %d\n",*_sta,*_sta_idx,iter_idx[*_sta_idx]);
    return read_info;
}

READ_t *bseq_read_map(bseq_file_t *fp, uint32_t batch, uint32_t *_n, uint32_t *_sta, uint32_t *_sta_idx, uint32_t *iter_map)
{
    READ_t *read_info = NULL;
    uint32_t idx = 0, base = 0, m = 0, seqi = *_sta, seqidx = *_sta_idx;
    int64_t kr1 = 1;
    kseq_t *seq1 = fp->ks;

    for (; (base < batch) && ((kr1 = kseq_read(seq1)) > 0); seqi++)
    {
        if (seqi == iter_map[seqidx])
        {
            if (idx >= m)
            {
                m = m == 0 ? 256 : m << 1;
                read_info = (READ_t *)realloc(read_info, m * sizeof(READ_t));
                if (read_info == NULL)
                {
                    fprintf(stderr, "[%s] calloc %ldGB read_info memory error!\n", __func__, m * sizeof(READ_t) / 1024 /1024 / 1024);
                    exit(1);
                }
            }
            read_info[idx].rid = seqi;
            read_info[idx].read_name = strdup(seq1->name.s);
            read_info[idx].read_length = seq1->seq.l;
            read_info[idx].read_seq = strdup(seq1->seq.s);
            base += read_info[idx++].read_length;
            seqidx++;
        }
    }

    *_n = idx;
    *_sta = seqi;
    *_sta_idx = seqidx;
    return read_info;
}