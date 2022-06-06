import os
import sys
import psutil
import numpy as np
from class_read import Maf
from class_read import Overlap
from class_read import Read_Info
from ove_in_maf_file import true_ove_in_maf_file

def load_maf(fp_maf, maf_list):
	with open(fp_maf, 'r', errors = 'ignore') as f:
		read_id = 0
		line = f.readline()#1
		while line:
			line = f.readline()#2
			ref = line.split()
			line = f.readline()#3
			read = line.split()
			maf = Maf(ref, read, read_id)
			if (read_id % 1000 == 0):
				print("load maf format of read ", read_id, end='\r')
			maf_list.append(maf)
			read_id = read_id + 1
			line = f.readline()#4
			line = f.readline()#0
	return read_id

def is_redundant_ove(ove_list, ove_n, dest_id, source_id):
	for k in range(0, ove_n[dest_id]):
		if (ove_list[dest_id][k].get_tid() == source_id):
			return 1
	return 0

# all reported ove
def load_ove(fp_ove, ove_list):
	with open(fp_ove,'r') as f:
		ove_num = 0
		line = f.readline()
		while line:
			line = line.split("\t")
			qid  = int(line[0])
			qlen = int(line[1])
			qsta = int(line[2])
			qend = int(line[3])
			rev = 0 if line[4] == '+' else 1
			tid  = int(line[5])
			tlen = int(line[6])
			tsta = int(line[7]) if rev == 0 else int(line[8])
			tend = int(line[8]) if rev == 0 else int(line[7])
			oveq = Overlap(qid, qlen, qsta, qend, rev, tid, tlen, tsta, tend)
			ove_list[qid].append(oveq)
			# print(qid, qlen, qsta, qend, rev, tid, tlen, tsta, tend)
			# print(tid, tlen, tsta, tend, rev, qid, qlen, qsta, qend)
			ove_num  = ove_num + 1
			if (qid % 1000 == 0):
				print("load ove of read ", qid, end='\r')
			line = f.readline()

	return ove_num

def count_true_positive_num(maf_list, ove_list, cor_read_all, wrg_read_all):
	true_positive_num = 0
	for i in range(0, len(maf_list)):
		if (i % 1000 == 0):
			print("processed read ", i, end='\r')
		for j in range(0, len(ove_list[i])):
			t_read_id = ove_list[i][j].get_tid()
			t_read_chr = maf_list[t_read_id].get_rchr()
			t_ref_sta = maf_list[t_read_id].get_refsta()
			t_ref_len = maf_list[t_read_id].get_reflen()
			t_read_len = maf_list[t_read_id].get_readlen()
			t_read_str = maf_list[t_read_id].get_readstr()

			t_Eindel = t_ref_len * 0.3
			if (t_read_str == '+'):
				t_ove_sta = min(ove_list[i][j].get_tsta(),ove_list[i][j].get_tend()) + t_ref_sta
				t_ove_end = max(ove_list[i][j].get_tend(),ove_list[i][j].get_tsta()) + t_ref_sta
			else:
				t_ove_end = t_ref_sta + t_ref_len - min(ove_list[i][j].get_tsta(),ove_list[i][j].get_tend())
				t_ove_sta = t_ref_sta + t_ref_len - max(ove_list[i][j].get_tsta(),ove_list[i][j].get_tend())
			
			q_read_id = ove_list[i][j].get_qid()
			q_read_chr = maf_list[q_read_id].get_rchr()
			q_ref_sta = maf_list[q_read_id].get_refsta()
			q_ref_len = maf_list[q_read_id].get_reflen()
			q_read_len = maf_list[q_read_id].get_readlen()
			q_read_str = maf_list[q_read_id].get_readstr()

			q_Eindel = q_ref_len * 0.3
			if (q_read_str == '+'):
				q_ove_sta = ove_list[i][j].get_qsta() + q_ref_sta
				q_ove_end = ove_list[i][j].get_qend() + q_ref_sta
			else:
				q_ove_end = q_ref_sta + q_ref_len - ove_list[i][j].get_qsta()
				q_ove_sta = q_ref_sta + q_ref_len - ove_list[i][j].get_qend()

			if (q_read_chr == t_read_chr and abs(q_ove_sta - t_ove_sta) < 2*(q_Eindel + t_Eindel) and abs(q_ove_end - t_ove_end) < 2*(q_Eindel + t_Eindel)):
				true_positive_num = true_positive_num + 1
				cor_read_all[q_read_id] = 1
				cor_read_all[t_read_id] = 1
			else:
				wrg_read_all[q_read_id] = 1
				wrg_read_all[t_read_id] = 1
				# if (q_read_id == 3981):
				# print(q_read_id, q_read_str, ove_list[i][j].get_qlen(), ove_list[i][j].get_qsta(), ove_list[i][j].get_qend(),'\t',t_read_id,t_read_str,ove_list[i][j].get_tlen(),ove_list[i][j].get_tsta(),ove_list[i][j].get_tend())

	return true_positive_num

def estimate_accurcy(reported_ove_num, true_positive_num):
	ove_accurcy = (true_positive_num / reported_ove_num) * 100
	print("overlap accurcy: " + "%0.5f" % ove_accurcy + "%(" + "%d" %true_positive_num + "/" + "%d" %reported_ove_num + ")")

def estimate_read_performance(maf_list, cor_read_all, wrg_read_all, read_all_num, read_N_num):
	read_all_num = read_all_num - read_N_num

	cor_read_all_num = 0
	wrg_read_part_num = 0
	wrg_read_all_num = 0
	nof_read_all_num = 0
	for i in range(0,len(cor_read_all)):
		# if (i == 985):
		# 	print('read %d' %i,'S%d' %(maf_list[i].get_rchr()+1)+'_%d'%(maf_list[i].get_rid()+1))
		# 	for k in range(0, len(maf_list[i])):
		# 		ref_sta = maf_list[i][k].get_refsta()
		# 		ref_end = maf_list[i][k].get_refsta() + maf_list[i][k].get_reflen()
		# 		if (maf_list[i].get_rchr() == maf_list[i][k].get_rchr() and ref_sta < maf_list[i].get_refsta() + maf_list[i].get_reflen() and ref_end > maf_list[i].get_refsta()):
		# 			print(k,'S%d' %(maf_list[i][k].get_rchr()+1)+'_%d'%(maf_list[i][k].get_rid()+1))
		if (cor_read_all[i] == 1 and wrg_read_all[i] == 0):
			cor_read_all_num = cor_read_all_num + 1
		elif (cor_read_all[i] == 1 and wrg_read_all[i] == 1):
			wrg_read_part_num = wrg_read_part_num + 1
		elif (cor_read_all[i] == 0 and wrg_read_all[i] == 1):
			wrg_read_all_num = wrg_read_all_num + 1
			# print(i,'S%d' %(maf_list[i].get_rchr()+1)+'_%d'%(maf_list[i].get_rid()+1), maf_list[i].get_readlen())
		else:
			nof_read_all_num = nof_read_all_num + 1
			# print('read %d' %i,'S%d' %(maf_list[i].get_rchr()+1)+'_%d'%(maf_list[i].get_rid()+1))
			# for k in range(0, len(maf_list[i])):
			# 	ref_sta = maf_list[i][k].get_refsta()
			# 	ref_end = maf_list[i][k].get_refsta() + maf_list[i][k].get_reflen()
			# 	if (maf_list[i].get_rchr() == maf_list[i][k].get_rchr() and ref_sta < maf_list[i].get_refsta() + maf_list[i].get_reflen() and ref_end > maf_list[i].get_refsta()):
			# 		print(k,'S%d' %(maf_list[i][k].get_rchr()+1)+'_%d'%(maf_list[i][k].get_rid()+1))
    					
	read_recall_all = (cor_read_all_num / read_all_num) * 100
	print("read recall(with TP overlaps): %0.5f" %read_recall_all + "%(" + "%d" %cor_read_all_num + "/" + "%d" %read_all_num + ") reads have found all correct overlaps")

	read_recall_part_wrg = (wrg_read_part_num / read_all_num) * 100
	print("read recall(with wrong overlaps): %0.5f" %read_recall_part_wrg + "%(" + "%d" %wrg_read_part_num + "/" + "%d" %read_all_num + ") reads have found part wrong overlaps")

	read_recall_wrg = (wrg_read_all_num / read_all_num) * 100
	print("read recall(with wrong overlaps): %0.5f" %read_recall_wrg + "%(" + "%d" %wrg_read_all_num + "/" + "%d" %read_all_num + ") reads have found all wrong overlaps")

	nof_read_all_num = nof_read_all_num - read_N_num
	read_recall_nof = (nof_read_all_num / read_all_num) * 100
	print("read recall(cannot find overlaps): %0.5f" %read_recall_nof + "%(" + "%d" %nof_read_all_num + "/" + "%d" %read_all_num + ") reads haven't found overlaps")

def estimate_transitive(fp_true_ove, fp_trans_ove, ove_list, true_ove_list, read_all_num):
	with open (fp_true_ove, 'r') as fp:
		line = fp.readline()
		while line:
			line = line.split("\t")
			qid = int(line[0])
			qlen = int(line[1])
			qsta = int(line[2])
			qend = int(line[3])
			rev  = int(line[4])
			tid = int(line[5])
			tlen = int(line[6])
			tsta = int(line[7])
			tend = int(line[8][:-1]) # get rid of \n
			ove  = Overlap(qid, qlen, qsta, qend, rev, tid, tlen, tsta, tend)
			true_ove_list[qid].append(ove)
			if (qid % 1000 == 0):
				print("load true ove of read ", qid, end='\r')
			line = fp.readline()
	
	with open(fp_trans_ove, 'w') as f:
		ove_n = [len(ove_list[i]) for i in range (read_all_num)]
		trans_ove_num = 0
		true_ove_num = 0
		same_ove = 0
		lack_ove = 0
		extra_ove = 0

		ove_n = [len(ove_list[i]) for i in range (read_all_num)]
		for i in range(0, read_all_num):
			for j in range(0, ove_n[i]):
				ovet = Overlap(ove_list[i][j].get_tid(), ove_list[i][j].get_tlen(), ove_list[i][j].get_tsta(), ove_list[i][j].get_tend(), ove_list[i][j].get_rev(), ove_list[i][j].get_qid(), ove_list[i][j].get_qlen(), ove_list[i][j].get_qsta(), ove_list[i][j].get_qend())
				if (is_redundant_ove(ove_list, ove_n, ove_list[i][j].get_tid(), ove_list[i][j].get_qid()) != 1):
					ove_list[ove_list[i][j].get_tid()].append(ovet)
		ove_n = [len(ove_list[i]) for i in range (read_all_num)]

		for i in range(0, read_all_num):
			if (i % 1000 == 0):
				print("processed read ", i, end='\r')
			if (len(ove_list[i]) != 0):
				visited = [0 for ii in range (read_all_num)]
				stack = []
				read_id = i
				read_len = ove_list[i][0].get_qlen()
				cov_q_sta_pos = 0
				cov_q_end_pos = 0
				last_read = Read_Info(read_id, 0, 0, 0, 0, 0)
				stack.append(last_read)
				visited[last_read.get_rid()] = 1
				while(len(stack) > 0):
					last_read = stack.pop()
					for j in range(0, ove_n[last_read.get_rid()]):
						if (visited[ove_list[last_read.get_rid()][j].get_tid()] == 1):
							continue
						if (last_read.get_rev() == 1):
							cov_q_sta_pos = last_read.get_sta_pos() - (ove_list[last_read.get_rid()][j].get_qend() - last_read.get_te())
							cov_q_end_pos = last_read.get_end_pos() - (ove_list[last_read.get_rid()][j].get_qsta() - last_read.get_ts())
						else:
							cov_q_sta_pos = last_read.get_sta_pos() + ove_list[last_read.get_rid()][j].get_qsta() - last_read.get_ts()
							cov_q_end_pos = last_read.get_end_pos() + ove_list[last_read.get_rid()][j].get_qend() - last_read.get_te()
						if (cov_q_sta_pos > ove_list[read_id][0].get_qlen() or cov_q_end_pos < 0):
							continue
						# if (read_id == 3981):
						# 	print(ove_list[last_read.get_rid()][j].get_qid(), ove_list[last_read.get_rid()][j].get_qlen(), ove_list[last_read.get_rid()][j].get_qsta(), ove_list[last_read.get_rid()][j].get_qend(), ove_list[last_read.get_rid()][j].get_rev(), ove_list[last_read.get_rid()][j].get_tid(), ove_list[last_read.get_rid()][j].get_tlen(), ove_list[last_read.get_rid()][j].get_tsta(), ove_list[last_read.get_rid()][j].get_tend())
						# 	print(last_read.get_rid(), last_read.get_sta_pos(), last_read.get_end_pos(), last_read.get_ts(), last_read.get_te(), last_read.get_rev())
						# 	print(cov_q_sta_pos, cov_q_end_pos)
						if (last_read.get_rid() != read_id):
							ove_trans = Overlap(read_id, read_len, cov_q_sta_pos, cov_q_end_pos, last_read.get_rev()^ove_list[last_read.get_rid()][j].get_rev(), ove_list[last_read.get_rid()][j].get_tid(), ove_list[last_read.get_rid()][j].get_tlen(), ove_list[last_read.get_rid()][j].get_tsta(), ove_list[last_read.get_rid()][j].get_tend())
							ove_list[read_id].append(ove_trans)
						cur_read = Read_Info(ove_list[last_read.get_rid()][j].get_tid(), cov_q_sta_pos, cov_q_end_pos, ove_list[last_read.get_rid()][j].get_tsta(), ove_list[last_read.get_rid()][j].get_tend(), last_read.get_rev()^ove_list[last_read.get_rid()][j].get_rev())
						stack.append(cur_read)
						visited[cur_read.get_rid()] = 1
				# if (len(ove_list[read_id]) > 1000):
				# if (read_id == 7000):
					# break

				# f.write('%d\t' %i + '%d\n' %len(ove_list[i]))
				# for j in range(0, len(ove_list[i])):
				# 	f.write('%d\t' %ove_list[i][j].get_qid() + '%d\t' %ove_list[i][j].get_qlen() + '%d\t' %ove_list[i][j].get_qsta() + '%d\t' %ove_list[i][j].get_qend() + '%d\t' %ove_list[i][j].get_rev() + '%d\t' %ove_list[i][j].get_tid() + '%d\t' %ove_list[i][j].get_tlen() + '%d\t' %ove_list[i][j].get_tsta() + '%d\n' %ove_list[i][j].get_tend())

			trans_ove_num += len(ove_list[i])
			edge_list = set()
			for j in range (0, len(ove_list[i])):
				edge_list.add(ove_list[i][j].get_tid())

			true_ove_num += len(true_ove_list[i])
			true_edge_list = set()
			for j in range (0, len(true_ove_list[i])):
				true_edge_list.add(true_ove_list[i][j].get_tid())

			edges = sorted(list(edge_list))
			true_edges = sorted(list(true_edge_list))
			# print(edges)
			# print(true_edges)
			# j for edges's idx, k for true_edges's idx
			j = 0
			k = 0
			while(k < len(true_edges) and j < len(edges)):
				if (edges[j] == true_edges[k]):
					# print('same', edges[j]) # intersection of edges and ture_edges
					j = j + 1
					k = k + 1
					same_ove = same_ove + 1
				elif (edges[j] > true_edges[k]): # edges lack of tid
					# print('lack', true_edges[k])
					k = k + 1
					lack_ove = lack_ove + 1
				elif (edges[j] < true_edges[k]): # edges has one extra tid
					# print('extra', edges[j])
					j = j + 1
					extra_ove = extra_ove + 1
			# if (i == 10):
			# 	break
			if (k < len(true_edges)):
				lack_ove += (len(true_edges) - k)
			if (j < len(edges)):
				extra_ove += (len(edges) - j)

			del ove_list[i][ove_n[i]:]

		print('true overlap num: %d; ' %true_ove_num + 'totally find %d overlaps(reported and transtive)' %trans_ove_num)
		print('overlaps true: %d' %same_ove)
		print('overlaps not find: %d' %lack_ove)
		print('overlaps extra(wrong): %d' %extra_ove)

def main():
	fp_ove = sys.argv[1]
	fp_maf = sys.argv[2]
	read_N_num = int(sys.argv[3])
	print('[Ove]evaluating dataset ' + fp_maf.split('/')[5].split('.')[0])
	fp_trans_ove = fp_maf.split('/')[5].split('.')[0] + '.trans_oves.paf'
	fp_true_ove = fp_maf.split('/')[5].split('.')[0] + '.true_oves.paf'

	maf_list = []
	read_all_num = load_maf(fp_maf, maf_list)
	print('memory usage: %.4f GB' %(psutil.Process(os.getpid()).memory_info().rss/1024/1024/1024))
	ove_list = [[] for i in range (read_all_num)]
	reported_ove_num = load_ove(fp_ove, ove_list)
	print('memory usage: %.4f GB' %(psutil.Process(os.getpid()).memory_info().rss/1024/1024/1024))

	cor_read_all = [0 for i in range (read_all_num)]
	wrg_read_all = [0 for i in range (read_all_num)]
	true_positive_num = count_true_positive_num(maf_list, ove_list, cor_read_all, wrg_read_all)
	estimate_accurcy(reported_ove_num, true_positive_num)
	estimate_read_performance(maf_list, cor_read_all, wrg_read_all, read_all_num, read_N_num)
	print('memory usage: %.4f GB' %(psutil.Process(os.getpid()).memory_info().rss/1024/1024/1024))

	# true_ove_list = [[] for i in range (read_all_num)]
	# estimate_transitive(fp_true_ove, fp_trans_ove, ove_list, true_ove_list, read_all_num)
	# print('memory usage: %.4f GB' %(psutil.Process(os.getpid()).memory_info().rss/1024/1024/1024))

	# cor_read_all = [0 for i in range (read_all_num)]
	# wrg_read_all = [0 for i in range (read_all_num)]
	# true_positive_num = count_true_positive_num(maf_list, ove_list, cor_read_all, wrg_read_all)
	# estimate_accurcy(reported_ove_num, true_positive_num)
	# estimate_read_performance(maf_list, cor_read_all, wrg_read_all, read_all_num, read_N_num)

if __name__ == "__main__":
	main()
