#
# This script generates a database of wormbase ids.
#


import gzip
import os.path
from peewee import *
from subprocess import check_output
from collections import OrderedDict
import re
from pprint import pprint as pp


build = "WS245"

gff_url = "ftp://ftp.wormbase.org/pub/wormbase/releases/{build}/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.{build}.annotations.gff3.gz".format(build = build)
gff = "c_elegans.{build}.gff3.gz".format(build = build)

if not os.path.exists(gff):
    comm = "curl {gff_url} > {gff}".format(**locals())
    print(check_output(comm, shell = True))


try:
    os.remove("celegans_gff.db")
except:
    pass




db = SqliteDatabase('celegans_gff.db')
db.connect() 

class BaseModel(Model):
    class Meta:
        database = db

class Feature(BaseModel):
  chrom = CharField()
  source = CharField()
  type_of = CharField()
  start = IntegerField()
  end = IntegerField()
  score = IntegerField(null = True)
  strand = FixedCharField(null = True, max_length = 1)
  reading_frame = IntegerField(null = True)

  class Meta:
    order_by = ('chrom',)

class Attribute(BaseModel):
  feature = ForeignKeyField(Feature, related_name="ID", null=True )
  id_key = CharField(index = True, null=True)
  id_value = CharField(index = True, null=True)
  
  class Meta:
    database = db

db.create_tables([Feature,Attribute])
db.execute_sql("CREATE VIEW IF NOT EXISTS feature_id AS SELECT chrom, source, type_of, start, end, score, strand, reading_frame, id_key, id_value FROM feature JOIN attribute ON feature.id == attribute.feature_id;")


def boolify(s):
    if s == 'True':
        return True
    if s == 'False':
        return False
    raise ValueError("huh?")

def autoconvert(s):
    for fn in (boolify, int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    return s

correct_ids = ["locus", "Name", "sequence_name"]
header = ["chrom", "source", "type_of", "start", "end", "score", "strand", "reading_frame", "attributes"]

with gzip.open(gff, 'rb') as f:
    c = 0
    while True:
        with db.atomic():
            try:
                line = f.next().strip().split("\t")
            except:
                break
            if line[0].startswith("#"):
                continue
            if line[1] == "WormBase":
                c += 1
                if c % 5000 == 0:
                    print c, record["chrom"], record["start"]
                #Put in bulk, every 1000th iteration
                record = OrderedDict(zip(header, map(autoconvert, line)))
                attributes = dict([(e.split("=")[0], autoconvert(e.split("=")[1])) for e in record["attributes"].split(";")])
                del record["attributes"]
                if record["score"] == '.':
                    record["score"] = None
                if record["reading_frame"] == ".":
                    record["reading_frame"] = None
                feature_id, success = Feature.get_or_create(**record)
                if success:
                    for k,v in attributes.items():
                        Attribute.get_or_create(feature = feature_id.get_id(), id_key = k, id_value = v)

db.close()

