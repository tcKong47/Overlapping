#load .maf file...
class Maf(object):
	def __init__(self, ref, read, rid):
		self.__refsta  = int(ref[2])
		self.__refend  = int(ref[2]) + int(ref[3])
		self.__reflen  = int(ref[3])
		self.__refstr  = ref[4]
		self.__rid     = rid
		self.__rchrid  = int(read[1].split('_')[1])-1
		self.__rchr    = int(read[1].split('_')[0][1:])-1
		self.__readlen = int(read[3])
		self.__readstr = read[4]

	def get_refsta(self):
		return self.__refsta

	def get_refend(self):
		return self.__refend

	def get_reflen(self):
		return self.__reflen

	def get_refstr(self):
		return self.__refstr

	def get_rid(self):
		return self.__rid

	def get_rchrid(self):
		return self.__rchrid

	def get_rchr(self):
		return self.__rchr

	def get_readlen(self):
		return self.__readlen

	def get_readstr(self):
		return self.__readstr

#load overlaps...
class Overlap(object):
	def __init__(self, qid, qlen, qsta, qend, rev, tid, tlen, tsta, tend):
		self.__qid  = qid
		self.__qlen = qlen
		self.__qsta = qsta
		self.__qend = qend
		self.__rev  = rev

		self.__tid  = tid
		self.__tlen = tlen
		self.__tsta = tsta
		self.__tend = tend

		# self.__mbp = int(ove[9])
		# self.__mln = int(ove[10])
		# self.__score = float(ove[11])

	def get_qid(self):
		return self.__qid

	def get_qlen(self):
		return self.__qlen

	def get_qsta(self):
		return self.__qsta

	def get_qend(self):
		return self.__qend

	def get_rev(self):
		return self.__rev

	def get_tid(self):
		return self.__tid

	def get_tlen(self):
		return self.__tlen

	def get_tsta(self):
		return self.__tsta

	def get_tend(self):
		return self.__tend
	
	# def get_mbp(self):
	# 	return self.__mbp

	# def get_mln(self):
	# 	return self.__mln

	# def get_score(self):
	# 	return self.__score


#load overlaps...
class Read_Info(object):
	def __init__(self, rid, sta_pos, end_pos, ts, te, rev):
		self.__rid  = rid
		self.__sta_pos = sta_pos
		self.__end_pos = end_pos
		self.__ts = ts
		self.__te  = te
		self.__rev  = rev

	def get_rid(self):
		return self.__rid

	def get_sta_pos(self):
		return self.__sta_pos

	def get_end_pos(self):
		return self.__end_pos

	def get_ts(self):
		return self.__ts

	def get_te(self):
		return self.__te

	def get_rev(self):
		return self.__rev