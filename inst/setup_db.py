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
    os.remove("elegans_gff.db")
except:
    pass


db = SqliteDatabase('celegans_gff.db')
db.connect() 

class BaseModel(Model):
    class Meta:
        database = db

class Region(BaseModel):
  chrom= CharField()
  start = IntegerField()
  end = IntegerField()
  #strand = CharField()
  #reading_frame = CharField()

  class Meta:
    order_by = ('chrom',)

class Id_Set(BaseModel):
  region = ForeignKeyField(Region, related_name="ID", null=True )
  type_of = CharField(index = True, null = True)
  id_key = CharField(index = True, null=True)
  id_value = CharField(index = True, null=True)
  
  class Meta:
    database = db
    indexes = (
            # create a unique on from/to/date
            (('region', 'type_of', 'id_key', 'id_value'), True),
            )

db.create_tables([Region,Id_Set])
db.execute_sql("CREATE VIEW IF NOT EXISTS region_id AS SELECT chrom, start, end, type_of, id_key, id_value FROM region JOIN id_set ON region.id == id_set.region_id;")


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
                if any([record["attributes"].find(k) > 0 for k in correct_ids]):
                    region_save = {k:v for k,v in record.items() if k in ["chrom","start","end"]}
                    region_id, success = Region.get_or_create(**region_save)
                    if not success:
                        region_id.save()
                    id_group = []
                    id_record = {}
                    for ID, VAL in re.findall("([^=;]+)=([^; ]+)", record["attributes"]):
                        if ID in correct_ids:
                            if VAL.find(":") > 0:
                                VAL = VAL.split(":")[1]
                            id_record["id_key"]= ID
                            id_record["type_of"]= record["type_of"]
                            id_record["id_value"] = VAL
                            id_record["region"] = region_id.get_id()
                            rec, success = Id_Set.get_or_create(**id_record)
                            rec.save()

db.close()

