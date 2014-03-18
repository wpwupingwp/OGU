from __future__ import unicode_literals
from django.db import models

class Main(models.Model):
    Taxon = models.IntegerField(db_column='Taxon') 
    Organism = models.TextField(db_column='Organism') 
    Accession = models.TextField(db_column='Accession') 
    Name = models.TextField(db_column='Name') 
    Type = models.TextField(db_column='Type') 
    Head = models.IntegerField(db_column='Head') 
    Tail = models.IntegerField(db_column='Tail') 
    Strand = models.TextField(db_column='Strand') 
    Sequence = models.TextField(db_column='Sequence') 
    Date = models.TextField(db_column='Date') 
    ID = models.IntegerField(db_column='ID',primary_key=True)
    class Meta:
        db_table = 'main'
        ordering=['ID']

