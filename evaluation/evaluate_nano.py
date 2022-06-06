import sys
import numpy as np

ove_list = []
read_list = []

class Overlap(object):
	def __init__(self, ove):
		self.__qchr = 0
		self.__qid  = int(ove[0])
		self.__qlen = int(ove[1])
		self.__qsta = int(ove[2])
		self.__qend = int(ove[3])
		self.__qstr = ove[4]

		self.__tchr = 0
		self.__tid  = int(ove[5])
		self.__tlen = int(ove[6])
		self.__tsta = int(ove[7])
		self.__tend = int(ove[8])
		
		self.__mbp = int(ove[9])
		self.__mln = int(ove[10])
		self.__score = float(ove[11])

	def get_qid(self):
    		return self.__qid

	def get_qchr(self):
		return self.__qchr

	def get_qlen(self):
		return self.__qlen

	def get_qsta(self):
		return self.__qsta

	def get_qend(self):
		return self.__qend

	def get_qstr(self):
		return self.__qstr

	def get_tid(self):
		return self.__tid

	def get_tchr(self):
		return self.__tchr

	def get_tlen(self):
		return self.__tlen

	def get_tsta(self):
		return self.__tsta

	def get_tend(self):
		return self.__tend
	
	def get_mbp(self):
		return self.__mbp

	def get_mln(self):
		return self.__mln

	def get_score(self):
		return self.__score

# paf file
def load_ove(fp_ove, ove_list):
	with open(fp_ove,'r') as f:
		line = f.readline() # line 0
		while line:
			line = line.split("\t")
			ove = Overlap(line)
			ove_list.append(ove)
			line = f.readline() # line 0

class Read(object):
	def __init__(self, ref, pos, rid, strand, head_bp, aligned_bp, tail_bp):
		self.__chr = ref
		self.__pos = pos
		self.__rid = rid
		self.__strand  = strand
		self.__readlen = head_bp + aligned_bp + tail_bp

	def get_chr(self):
		return self.__chr

	def get_pos(self):
		return self.__pos

	def get_rid(self):
		return self.__rid

	def get_strand(self):
		return self.__strand

	def get_readlen(self):
		return self.__readlen

# fasta/fastq file
def load_fa(fp_fa, read_list):
	with open(fp_fa, 'r') as f:
		read_id = 0
		line = f.readline() # line 0
		while line:
			if (line[0] == '>'): # fasta file
				line = line.split('_') # >NC-001148_239450;aligned_0_R_20_1212_49
				ref = line[0][1:] # get rid of '>'
				pos = int(line[1].split(';')[0]) # get rid of 'aligned' or 'unaligned'
				rid = read_id # rid not use index int(line[2]) in file
				strand = line[3]
				head_bp = int(line[4]) # convert string to int
				aligned_bp = int(line[5]) # convert string to int
				tail_bp = int(line[6][:-1]) # get rid of '\n'
				read = Read(ref, pos, rid, strand, head_bp, aligned_bp, tail_bp)
				read_list.append(read)
				line = f.readline() # line 1
				read_id = read_id + 1
				line = f.readline() # line 0
				# break
			elif (line[0] == '@'): # fastq file
				pass

def main():
	fp_ove = sys.argv[1]
	fp_fa = sys.argv[2]
	load_ove(fp_ove, ove_list)
	load_fa(fp_fa, read_list)
	read_N_num = int(sys.argv[3])

	wrg_read_all = [0 for i in range (len(read_list))]
	cor_read_all = [0 for i in range (len(read_list))]

	true_positive_num = count_true_positive_num(read_list, ove_list, wrg_read_all, cor_read_all)
	reported_ove_num = len(ove_list)
	estimate_accurcy(reported_ove_num, true_positive_num)
	estimate_read_performance(wrg_read_all, cor_read_all, read_N_num)


def count_true_positive_num(read_list, ove_list, wrg_read_all, cor_read_all):
	true_positive_num = 0
	for i in range(0, len(ove_list)):
		t_read_id = ove_list[i].get_tid()
		t_read_chr = read_list[t_read_id].get_chr()
		t_ref_sta = read_list[t_read_id].get_pos()
		t_read_len = read_list[t_read_id].get_readlen()
		t_read_str = read_list[t_read_id].get_strand()

		t_Eindel = t_read_len * 0.25
		if (t_read_str == 'F'):
			t_ove_sta = min(ove_list[i].get_tsta(), ove_list[i].get_tend()) + t_ref_sta
			t_ove_end = max(ove_list[i].get_tend(), ove_list[i].get_tsta()) + t_ref_sta
		else:
			t_ove_end = t_ref_sta + t_read_len - \
				min(ove_list[i].get_tsta(), ove_list[i].get_tend())
			t_ove_sta = t_ref_sta + t_read_len - \
				max(ove_list[i].get_tsta(), ove_list[i].get_tend())

		q_read_id = ove_list[i].get_qid()
		q_read_chr = read_list[q_read_id].get_chr()
		q_ref_sta = read_list[q_read_id].get_pos()
		q_read_len = read_list[q_read_id].get_readlen()
		q_read_str = read_list[q_read_id].get_strand()

		q_Eindel = q_read_len * 0.25
		if (q_read_str == 'F'):
			q_ove_sta = ove_list[i].get_qsta() + q_ref_sta
			q_ove_end = ove_list[i].get_qend() + q_ref_sta
		else:
			q_ove_end = q_ref_sta + q_read_len - ove_list[i].get_qsta()
			q_ove_sta = q_ref_sta + q_read_len - ove_list[i].get_qend()

		if (q_read_chr == t_read_chr and abs(q_ove_sta - t_ove_sta) < 2*(q_Eindel + t_Eindel) and abs(q_ove_end - t_ove_end) < 2*(q_Eindel + t_Eindel)):
			true_positive_num = true_positive_num + 1
			cor_read_all[q_read_id] = 1
			cor_read_all[t_read_id] = 1
		else:
			wrg_read_all[q_read_id] = 1
			wrg_read_all[t_read_id] = 1

	return true_positive_num

def estimate_accurcy(reported_ove_num, true_positive_num):
	ove_accurcy = (true_positive_num / reported_ove_num) * 100
	print("overlap accurcy: " + "%0.5f" % ove_accurcy + "%(" + "%d" %true_positive_num + "/" + "%d" %reported_ove_num + ")")

def estimate_read_performance(wrg_read_all, cor_read_all, read_N_num):
	read_all_num = len(read_list)
	read_all_num = read_all_num - read_N_num

	cor_read_all_num = 0
	wrg_read_part_num = 0
	wrg_read_all_num = 0
	nof_read_all_num = 0

	for i in range(0,len(cor_read_all)):
		if (cor_read_all[i] == 1 and wrg_read_all[i] == 0):
			cor_read_all_num = cor_read_all_num + 1
		elif (cor_read_all[i] == 1 and wrg_read_all[i] == 1):
			wrg_read_part_num = wrg_read_part_num + 1
		elif (cor_read_all[i] == 0 and wrg_read_all[i] == 1):
			wrg_read_all_num = wrg_read_all_num + 1
		else:
			nof_read_all_num = nof_read_all_num + 1

	read_recall_all = (cor_read_all_num / read_all_num) * 100
	print("read recall(with TP overlaps): %0.5f" %read_recall_all + "%(" + "%d" %cor_read_all_num + "/" + "%d" %read_all_num + ") reads have found all correct overlaps")

	read_recall_part_wrg = (wrg_read_part_num / read_all_num) * 100
	print("read recall(with wrong overlaps): %0.5f" %read_recall_part_wrg + "%(" + "%d" %wrg_read_part_num + "/" + "%d" %read_all_num + ") reads have found part wrong overlaps")

	read_recall_wrg = (wrg_read_all_num / read_all_num) * 100
	print("read recall(with wrong overlaps): %0.5f" %read_recall_wrg + "%(" + "%d" %wrg_read_all_num + "/" + "%d" %read_all_num + ") reads have found all wrong overlaps")

	nof_read_all_num = nof_read_all_num - read_N_num
	read_recall_nof = (nof_read_all_num / read_all_num) * 100
	print("read recall(cannot find overlaps): %0.5f" %read_recall_nof + "%(" + "%d" %nof_read_all_num + "/" + "%d" %read_all_num + ") reads haven't found overlaps")

if __name__ == "__main__":
	main()
