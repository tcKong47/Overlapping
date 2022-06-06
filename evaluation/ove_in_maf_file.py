import sys
from class_read import Maf
from class_read import Overlap

def load_maf(fp_maf, maf_list):
	with open(fp_maf, 'r', errors = 'ignore') as f:
		read_id = 0
		tmp_chr_id = -1
		line = f.readline()#1
		while line:
			line = f.readline()#2
			ref = line.split()
			line = f.readline()#3
			read = line.split()
			# print(read[:-1])
			maf = Maf(ref, read, read_id)
			chr_id = int(read[1].split('_')[0][1:])-1
			if (chr_id != tmp_chr_id):
				maf_list.append([])
				tmp_chr_id = chr_id
			maf_list[chr_id].append(maf)
			read_id = read_id + 1
			line = f.readline()#4
			line = f.readline()#0
	return read_id

def processing_maf(fp_true_ove, maf_list, ove_list, read_all_num):
	with open(fp_true_ove, 'w') as f:
		# sorting maf_list using key refsta
		for i in range(0, len(maf_list)):
			# print(len(maf_list[i]))
			maf_list[i].sort(key = lambda x: x.get_refsta())
			for j in range(0, len(maf_list[i])):
				# print(i,j)
				refsta1  = maf_list[i][j].get_refsta()
				refend1  = maf_list[i][j].get_refend()
				reflen1  = maf_list[i][j].get_reflen()
				rid1     = maf_list[i][j].get_rid()
				readlen1 = maf_list[i][j].get_readlen()
				readstr1 = maf_list[i][j].get_readstr()
				for k in range(j + 1, len(maf_list[i])):
					refsta2  = maf_list[i][k].get_refsta()
					refend2  = maf_list[i][k].get_refend()
					reflen2  = maf_list[i][k].get_reflen()
					rid2     = maf_list[i][k].get_rid()
					readlen2 = maf_list[i][k].get_readlen()
					readstr2 = maf_list[i][k].get_readstr()
					if (refsta2 < refend1 - 500):
						qid  = rid1
						qlen = reflen1
						qsta = refsta1
						qend = refend1
						tid  = rid2
						tlen = reflen2
						tsta = refsta2
						tend = refend2
						rev  = 0 if (readstr1 == readstr2) else 1
						oveq = Overlap(qid, qlen, qsta, qend, rev, tid, tlen, tsta, tend)
						ove_list[qid].append(oveq)
						ovet = Overlap(tid, tlen, tsta, tend, rev, qid, qlen, qsta, qend)
						ove_list[tid].append(ovet)
						# print(qid, qlen, qsta, qend, rev, tid, tlen, tsta, tend)
						# print(tid, tlen, tsta, tend, rev, qid, qlen, qsta, qend)
					else:
						break

		for i in range(0, read_all_num):
			# f.write('%d\t' %i + '%d\n' %len(ove_list[i]))
			for j in range(0, len(ove_list[i])):
				f.write('%d\t' %ove_list[i][j].get_qid() + '%d\t' %ove_list[i][j].get_qlen() + '%d\t' %ove_list[i][j].get_qsta() + '%d\t' %ove_list[i][j].get_qend() + '%d\t' %ove_list[i][j].get_rev() + '%d\t' %ove_list[i][j].get_tid() + '%d\t' %ove_list[i][j].get_tlen() + '%d\t' %ove_list[i][j].get_tsta() + '%d\n' %ove_list[i][j].get_tend())

def true_ove_in_maf_file(fp_true_ove, fp_maf):
	maf_list = []
	read_all_num = load_maf(fp_maf, maf_list)
	ove_list = [[] for i in range (read_all_num)]
	processing_maf(fp_true_ove, maf_list, ove_list, read_all_num)
	return ove_list

if __name__ == "__main__":
	fp_maf = sys.argv[1]
	fp_true_ove = fp_maf.split('.')[0] + '.true_oves.paf'
	print(fp_true_ove)
	true_ove_in_maf_file(fp_true_ove, fp_maf)