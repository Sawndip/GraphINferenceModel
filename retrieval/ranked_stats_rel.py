#/usr/bin/env python
#
# Takes a TREC style results and qrel file and highlight which are the relevant docs

import fileinput
from qrel_reader import *

qrels = QrelReader("params/medtrack-all.qrel").get_qrels()

for line in fileinput.input():
	qId = line.split()[0]

	if qId not in qrels:
		print qId, "--"
	else:
		doc = line.split()[2]
		rel = "*" if doc in qrels[qId] else " "
		print "%s %s" % (line.strip(), rel)