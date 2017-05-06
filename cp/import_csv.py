#!/usr/bin/python3

import sqlite3


con = sqlite3.connect('data.sql')
cur = con.cursor()
cur.execute('''DROP TABLE IF EXISTS DNABankv3;''')
cur.execute('''CREATE TABLE DNABankv3 (BOP TEXT PRIMARY KEY,Type1 TEXT ,WeiZhi TEXT
,Type2 TEXT ,CanZha TEXT ,Chinese TEXT ,FamilyNumber TEXT ,Family TEXT ,Genus
TEXT ,Species TEXT ,Subspecies TEXT ,Country TEXT ,Province TEXT ,City TEXT
,Location TEXT ,Latitude TEXT ,Longitude TEXT ,Place TEXT ,Altitude TEXT
,Collector TEXT ,CollectNumber TEXT ,CollectDate TEXT ,Barcode TEXT ,Sampler
TEXT ,ChineseFamily TEXT ,ChineseSubFamily TEXT ,Info TEXT ,item1 TEXT ,item2
TEXT ,item3);''')
with open('./DNABank-v4.csv', 'r') as raw:
    for line in raw:
        line = line.strip().split(',')
        if len(line) != 30:
            line.extend(['']*(30-len(line)))
        command = '''INSERT INTO DNABankv3 (BOP, Type1, WeiZhi, Type2, CanZha,
        Chinese, FamilyNumber, Family, Genus, Species, Subspecies, Country,
        Province, City, Location, Latitude, Longitude, Place, Altitude,
        Collector, CollectNumber, CollectDate, Barcode, Sampler, ChineseFamily,
        ChineseSubFamily, Info, item1, item2, item3) VALUES (?,?,?,?,?,?,?,?,?,
        ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?, ?,?,?);'''
        try:
            cur.execute(command, line[:30])
        except:
            print(line)
cur.execute('''CREATE UNIQUE INDEX BOP_ID ON DNABankv3 (BOP);''')
cur.close()
# if miss, then every change miss
con.commit()
con.close()
