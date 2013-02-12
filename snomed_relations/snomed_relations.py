#!/usr/bin/env python
#
# snomed_relations
#
# Read the sct1_Relationships_Core_INT_20110131.txt file from stdin and return the relations for a given CUI
#


import sqlite3, fileinput, os, sys

BLACK_LIST = ['362955004', # Inactive concept (inactive concept)
			  '363660007', # Ambiguous concept (inactive concept)
			  '363661006', # Reason not stated concept (inactive concept)
			  '363662004', # Duplicate concept (inactive concept)
			  '363663009', # Outdated concept (inactive concept)
			  '363664003', # Erroneous concept (inactive concept)
			  '370126003', # Moved elsewhere (inactive concept)
			  '443559000'] # Limited status concept (inactive concept)
 
DATA_DIR = os.environ["HOME"] + "/phd/papers/network_based_retrieval_as_inference/snomed_relations"
DB_FILE = DATA_DIR+'/snomed_rel.db'

def get_relationships(cui):
  conn = sqlite3.connect(DB_FILE)
  c = conn.cursor()
  key = (cui,)
  c.execute("select cui2, reltype, relcharacteristic from crel where cui1=?", key)
  rel = [(str(cui2[0].strip()), str(cui2[1].strip()), str(cui2[2].strip())) for cui2 in c]
  return rel
    
def build_relation_db():
  #os.system("sqlite3 %s < umls_rel.sql" % DB_FILE)
  print "Building DB", DB_FILE
  
  conn = sqlite3.connect(DB_FILE)
  c = conn.cursor()
  
  schema = '''CREATE TABLE crel (
       cui1 varchar(20),
       cui2 varchar(20),
       reltype varchar(20),
       relcharacteristic int,
       PRIMARY KEY(cui1, cui2)
);'''

  
  c.execute(schema)
  print "Database creating, now populating..."

  for (count, line) in enumerate(fileinput.input()):
    if count == 0:
      continue

    # see the field description from 
    # http://www.ihtsdo.org/fileadmin/user_upload/doc/tig/index.html
    # note that the reltype are snomed concepts in themselves
    values = line.split('\t')
    cui1 = values[1].lower()
    reltype = values[2]
    cui2 = values[3].lower()
    relcharacteristic = int(values[4])

    if cui1 not in BLACK_LIST and cui2 not in BLACK_LIST:
    	c.execute("insert or ignore into crel (cui1, cui2, reltype, relcharacteristic) values ('%s', '%s', '%s', '%i')" % (cui1, cui2, reltype, relcharacteristic))
    
    if count % 100000 == 0:
      print `count/100000` + "00K..",
      sys.stdout.flush()
      
  conn.commit()
  c.close()
  
if __name__ == "__main__":
  if os.path.exists(DB_FILE):
    cui = sys.argv[1]
    for (c, r, t) in get_relationships(cui):
      print cui, c, r, t
  else:
    build_relation_db()