#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "ktime.h"
#include "main.h"
#include "layout.h"
#include "overlapping.h"

//global variable
int x_read;
int k;
int w;
int l;
float r_n;
int m_b;
int min_ove;
int max_hang;
float max_hang_ratio;
int top_n;
float cov_ratio;
int part_size;
int X_read;

int thread_n;
char *temp_file_perfix;
char temp_ove_dir[1024];
char temp_ove_link_dir[1024];
int waitingLen;
int batch_size_base;
int	batch_size_read;
double realtime0;

void init_map_param(param_map *opt)
{
	opt->x_read = 4;
	opt->k = 15;
	opt->w = 5;
	opt->l = 11;
	opt->r_n = 0.005;

	opt->m_b = 100;
	opt->min_ove = 500;
	opt->max_hang = 1500;
	opt->max_hang_ratio = 0.8;
	opt->top_n = 2;

	opt->cov_ratio = 0.33;
	opt->part_size = 3;
	opt->X_read = 10;

	opt->read_type = 1;
	opt->thread_n = 8;
	opt->batch_size = 100000;
}

int help_usage(param_map *opt)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:	overlapping\n");

	fprintf(stderr, "Usage:		ove <read.fa/fq> [options]\n");
	fprintf(stderr, "\nProgram options:   \n");
    fprintf(stderr, "    -h --help                           Help menu.\n");
    fprintf(stderr, "    -f <temporary file>        [STR]    The temporary file for storing overlapping information.\n");
	fprintf(stderr, "    -p --read-type             [INT]    Specifiy the type of reads and set multiple paramters unless overriden. [%d]\n", opt->read_type);
	fprintf(stderr, "                                        type 1: ont2d (Nanopore R10 reads): error rate 12%%\n");
	fprintf(stderr, "                                        type 2: ont new reads: error rate 1%%\n");
	fprintf(stderr, "    -t --thread-n              [INT]    thread n. [%d]\n", opt->thread_n);
	fprintf(stderr, "    -s --batch-size            [INT]    index batch size, query batch size *= 10. [%d]\n", opt->batch_size);

	fprintf(stderr, "\nAlgorithm options:\n");

	fprintf(stderr, "    -x --x-longest-read        [INT]    x percent longest reads of all reads. [%d]\n", opt->x_read);
	fprintf(stderr, "    -k --k-mer                 [INT]    k-mer length when indexing reads. [%d]\n", opt->k);
	fprintf(stderr, "    -w --window-size           [INT]    window size when indexing reads. [%d]\n", opt->w);
	fprintf(stderr, "    -l --l-mer                 [INT]    l-mer length for sorting minimizers. [%d]\n", opt->l);
	fprintf(stderr, "    -r --repeat-n              [FLOAT]  filter hits > repeat-n. [%.3f]\n", opt->r_n);

	fprintf(stderr, "    -b --matching-bases        [INT]    minimal matching bases. [%d]\n", opt->m_b);
	fprintf(stderr, "    -m --min-ove               [INT]    minimal overlap length. [%d]\n", opt->min_ove);
	fprintf(stderr, "    -a --max-hang              [INT]    maximal allowed overhang length. [%d]\n", opt->max_hang);
	fprintf(stderr, "    -T --max-hang-ratio        [FLOAT]  max overhang to mapping length ratio. [%.3f]\n", opt->max_hang_ratio);
	fprintf(stderr, "    -n --top-n                 [INT]    keeping top n overlaps. [%d]\n", opt->top_n);

	fprintf(stderr, "    -c --cov-ratio             [FLOAT]  coverage ratio for next iteration. [%.3f]\n", opt->cov_ratio);
	fprintf(stderr, "    -S --read-part-size        [INT]    size of read part whose cov is not enough for next iteration index. [%d]\n", opt->part_size);
	fprintf(stderr, "    -X --x-random-read         [INT]    X percent random reads with lacking coverage choosen for next iteration. [%d]\n", opt->X_read);
	return 1;
}

static const char *short_option = "f:p:t:s:x:k:w:l:r:b:m:a:T:n:c:S:X:h";

static struct option long_option[] = {
	{"temp-file-prefix", required_argument, NULL, 'f'},
	{"read-type", required_argument, NULL, 'p'},
	{"thread-n", required_argument, NULL, 't'},
	{"batch-size", required_argument, NULL, 's'},

	{"x-longest-read", required_argument, NULL, 'x'},
	{"k-mer", required_argument, NULL, 'k'},
	{"window-size", required_argument, NULL, 'w'},
	{"l-mer", required_argument, NULL, 'l'},
	{"repeat-n", required_argument, NULL, 'r'},

	{"matching-bases", required_argument, NULL, 'b'},
	{"min-ove", required_argument, NULL, 'm'},
	{"max-hang", required_argument, NULL, 'a'},
	{"max-hang-ratio", required_argument, NULL, 'T'},
	{"top-n", required_argument, NULL, 'n'},

	{"cov-ratio", required_argument, NULL, 'c'},
	{"read-part-size", required_argument, NULL, 'S'},
	{"x-random-read", required_argument, NULL, 'X'},
	{"help", no_argument, NULL, 'h'},
	{0, 0, 0, 0}};

int main(int argc, char *argv[])
{
	realtime0 = realtime();
	int c;
	float error, t, q, els;
	param_map *opt = (param_map* )calloc(1, sizeof(param_map));
	
	init_map_param(opt);

	while((c = getopt_long(argc, argv, short_option, long_option, NULL)) != -1)
	{
		switch(c)
		{
			case 'f': opt->temp_file_perfix = strdup(optarg); break;
			case 'p': opt->read_type = atoi(optarg); break;
			case 't': opt->thread_n = atoi(optarg); break;
			case 's': opt->batch_size = atoi(optarg); break;

			case 'x': opt->x_read = atoi(optarg); break;
			case 'k': opt->k = atoi(optarg); break;
			case 'w': opt->w = atoi(optarg); break;
			case 'l': opt->l = atoi(optarg); break;
			case 'r': opt->r_n = atof(optarg); break;

			case 'b': opt->m_b = atoi(optarg); break;
			case 'm': opt->min_ove = atoi(optarg); break;
			case 'a': opt->max_hang = atoi(optarg);break;
			case 'T': opt->max_hang_ratio = atof(optarg);break;
			case 'n': opt->top_n = atoi(optarg); break;

			case 'c': opt->cov_ratio = atof(optarg);break;
			case 'S': opt->part_size = atoi(optarg);break;
			case 'X': opt->X_read = atoi(optarg);break;
			case 'h': return help_usage(opt); break;
			default: return help_usage(opt); break;
		}
	} 

	if (argc - optind < 1)
		return help_usage(opt);

	if (opt->read_type > 2)
	{
		fprintf(stderr, "Input warning: read type should be 1 or 2, 1 default\n");
		opt->read_type = 1;
	}
	if (opt->thread_n < 1 || opt->thread_n > 32)
	{
		fprintf(stderr, "Input error: -t cannot be less than 1 or more than 32\n");
		exit(1);
	}
	// if (opt->batch_size < 100000 || opt->batch_size > 500000)
	// {
	// 	fprintf(stderr, "Input error: -s cannot be less than 1 or more than 32\n");
	// 	exit(1);
	// }
	if (opt->x_read <= 0 || opt->x_read > 100)
	{
		fprintf(stderr, "Input error: -x cannot be less than 0 or more than 100\n");
		exit(1);
	}
	if (opt->k < 10 || opt->k > 20)
	{
		fprintf(stderr, "Input error: -k cannot be less than 10 or more than 20\n");
		exit(1);
	}
	if (opt->w < 1 || opt->w > 10)
	{
		fprintf(stderr, "Input error: -w cannot be less than 1 or more than 10\n");
		exit(1);
	}
	if (opt->l < 6 || opt->l > 14)
	{
		fprintf(stderr, "Input error: -l cannot be less than 6 or more than 14\n");
		exit(1);
	}
	if (opt->l >= opt->k)
	{
		fprintf(stderr, "Input error: -l cannot be more than or equal to -k\n");
		exit(1);
	}
	if (opt->m_b < opt->k)
	{
		fprintf(stderr, "Input error: -b cannot be less than -k\n");
		exit(1);
	}
	if (opt->min_ove < opt->k)
	{
		fprintf(stderr, "Input error: -m cannot be less than -k\n");
		exit(1);
	}
	if (opt->max_hang < opt->k)
	{
		fprintf(stderr, "Input error: -a cannot be less than -k\n");
		exit(1);
	}
	if (opt->max_hang_ratio <= 0 || opt->max_hang_ratio > 1)
	{
		fprintf(stderr, "Input error: -T cannot be less than 0 or more than 1\n");
		exit(1);
	}
	if (opt->top_n < 1)
	{
		fprintf(stderr, "Input error: -n cannot be less than 1\n");
		exit(1);
	}
	if (opt->cov_ratio < 0.09 || opt->cov_ratio > 0.61)
	{
		fprintf(stderr, "Input error: -c cannot be less than 0.1 or more than 0.6\n");
		exit(1);
	}
	if (opt->part_size < 1 || opt->part_size > 10)
	{
		fprintf(stderr, "Input error: -S cannot be less than 1 or more than 10\n");
		exit(1);
	}
	if (opt->X_read < 1 || opt->X_read > 100)
	{
		fprintf(stderr, "Input error: -X cannot be less than 1 or more than 100\n");
		exit(1);
	}
	
	char *read_fastq;
	read_fastq = strdup(argv[optind]);
	char *index_fastq;
	index_fastq = strdup(argv[optind]);

	els = 0.05;
	error = 0.2;
	if (opt->read_type == 1) //ont 88%
		error = 0.12;
	else if (opt->read_type == 2) //ont 99%
		error = 0.01;

	t = log10(els)/log10(1 - pow(1 - error, opt->k));
	q = 1/error - opt->k * pow(1 - error, opt->k)/(1 - pow(1 - error, opt->k));
	waitingLen = (int)(t * q);

	x_read = opt->x_read;
	k = opt->k;
	w = opt->w;
	l = opt->l;
	r_n = opt->r_n;

	m_b = opt->m_b;
	min_ove = opt->min_ove;
	max_hang = opt->max_hang;
	max_hang_ratio = opt->max_hang_ratio;
	top_n = opt->top_n;

	cov_ratio = opt->cov_ratio;
	part_size = opt->part_size;
	X_read = opt->X_read;

	batch_size_read = 1000000;
	batch_size_base = 1000000000;
	thread_n = opt->thread_n;
	temp_file_perfix = opt->temp_file_perfix;

	memset(temp_ove_dir, 0, 1024);
	memset(temp_ove_link_dir, 0, 1024);
	if (temp_file_perfix == NULL)
	{
		strcpy(temp_ove_dir, "./oves.paf");
		strcpy(temp_ove_link_dir, "./oves_link.paf");
	}
	else
	{
		strcpy(temp_ove_dir, temp_file_perfix);
		strcat(temp_ove_link_dir, temp_file_perfix);
		strcat(temp_ove_dir, "_oves.paf");
		strcat(temp_ove_link_dir, "_oves_link.paf");
	}

	finding_overlapping(index_fastq, read_fastq, temp_ove_dir, temp_ove_link_dir);

	if (opt != NULL)	free(opt);
	if (read_fastq != NULL) {free(read_fastq); read_fastq = NULL;}
	if (index_fastq != NULL) {free(index_fastq); index_fastq = NULL;}

	fprintf(stderr, "[Result]\tReal time:%.3f sec\t%.3fh\tCPU:%.3f sec\tMemory peak:%.3f GB\n", realtime() - realtime0, (realtime() - realtime0) / 3600.0, cputime(), peak_memory() / 1024.0 / 1024.0);

	return 0;
}