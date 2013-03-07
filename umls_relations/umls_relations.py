#!/usr/bin/env python
#
# umls_relations
#
# Read the MRREL file from stdin and return the relations for a given CUI
#


import sqlite3, fileinput, os, sys

SOURCE = "ALL" # MTH = UMLS Meta, MSH = MeSh, SNOMEDCT
DATA_DIR = os.environ["HOME"] + "/phd/papers/network_based_retrieval_as_inference/umls_relations"
DB_FILE = DATA_DIR+'/umls_rel-%s.db' % SOURCE

def get_relationships(cui):
  conn = sqlite3.connect(DB_FILE)
  c = conn.cursor()
  key = (cui,)
  c.execute("select cui2, reltype from crel where cui1=?", key)
  rel = [(str(cui2[0].strip()), str(cui2[1].strip())) for cui2 in c]
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
       PRIMARY KEY(cui1, cui2)
);'''

  
  c.execute(schema)
  print "Database creating, now populating..."

  for (count, line) in enumerate(fileinput.input()):
    # see the field description from 
    # http://www.ncbi.nlm.nih.gov/books/NBK9685/table/ch03.T.related_concepts_file__mrrelrrf/?report=objectonly
    values = line.split('|')
    cui1 = values[0].lower()
    reltype = values[3]
    cui2 = values[4].lower()
    
    if (values[10] == SOURCE or SOURCE == "ALL"):
      c.execute("insert or ignore into crel (cui1, cui2, reltype) values ('%s', '%s', '%s')" % (cui1, cui2, reltype))
    
    if count % 1000000 == 0:
      print `count/1000000` + "M..",
      sys.stdout.flush()
      
  conn.commit()
  c.close()
  
if __name__ == "__main__":
  if os.path.exists(DB_FILE):
    cui = sys.argv[1]
    for (c, r) in get_relationships(cui):
      print cui, c, r
  else:
    build_relation_db()