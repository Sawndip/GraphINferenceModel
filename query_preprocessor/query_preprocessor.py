#!/usr/bin/env python
#
# Converts an Indri style query file into an extended SNOMED CT query

import argparse, commands, re
import odict

from iq_reader import *
from umls2snomed import *
from xml.sax.saxutils import escape

cui_regex = re.compile('([0-9]+)\s+(C[0-9]+):(.*) \[(.*)\]')
scui_desc_regex = re.compile('(.*)\((.*)\)')

mapper = UMLS2SNOMED()

def getNameAndSemanticType(scui):
	output = commands.getoutput('snomed_lookup %s' % scui)
	desc = output.split("\t")[2]
	(name, type) = scui_desc_regex.findall(desc)[0] 
	return (name.strip(), type)

def callMetaMap(txt, extendedFormat):
	cmd = 'echo %s | metamap11v2 -I -R SNOMEDCT -y' % txt
	output = commands.getoutput(cmd)
	
	phrase = ""
	phrases = odict.odict()
	senses = []
	cuis = []
	cands = True

	for (count, line) in enumerate(output.split('\n')):
		if line.startswith("Phrase"):
			if len(phrases) > 0:
				senses = []
				cuis = []
				cands = True
			phrase = line.split('"')[1]
			phrases[phrase] = senses
			del senses[0:len(senses)]
		if line.startswith("Meta Mapping"):
			score = re.findall('[0-9]+', line)[0]
			cands = False
			cuis = []
			senses.append((score, cuis))
			
		else:
			if not cands:
				cuis += [ (score, txt, type, cui) for (score, cui, txt, type) in cui_regex.findall(line) if cui not in [c[3] for c in cuis] ]
				
	if extendedFormat:
		for p in phrases:
			print "\t\t\t<phrase txt='%s'>" % p
			for (score, sense) in phrases[p]:
				print "\t\t\t\t<sense score='%s'>\n" % score,
				for concept in sense:
					# UMLS entry
					#print "\t\t\t\t\t<concept score='%s' name='%s' type='%s'>%s</concept>" % concept
					scuis = mapper.to_snomed(concept[3])
					if scuis != None:
						for scui in scuis.split():
							(name, type) = getNameAndSemanticType(scui)
							print "\t\t\t\t\t<concept score='%s' name=\"%s\" type='%s'>%s</concept>" % (concept[0], escape(name), escape(type), scui)
							#print scui,
					
				print "\t\t\t\t</sense>"
			print "\t\t\t</phrase>"
	else:
		print "\t\t\t",
		for p in phrases:
			for (score, sense) in phrases[p]:
				for concept in sense:
					
					scuis = mapper.to_snomed(concept[3])
					if scuis != None:
						for scui in scuis.split():
							print "%s_%s" % (scui, concept[0]),
		print ""

def handleQuery(id, txt, extendedFormat=False):
	callMetaMap(txt, extendedFormat)

if __name__ == "__main__":
	''' Main methods '''
	parser = argparse.ArgumentParser(description="Converts an Indri style query file into an extended SNOMED CT query")
	parser.add_argument("iq_file")
	parser.add_argument('-e', action='store_true', help="Print out in extended format")

	queries = IQReader(parser.parse_args().iq_file).get_queries()
	print '<parameters>'
	for (qid, txt) in queries.items():
		print '\t<query>'
		print '\t\t<number>%s</number>' % qid
		print '\t\t<text>'
		handleQuery(qid, txt.strip(), parser.parse_args().e)
		print '\t\t</text>'
		print '\t</query>'
	print '</parameters>'