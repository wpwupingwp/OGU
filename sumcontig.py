import re

for area in range(12):
    for cp in range(139):
        if cp<100:
            filename=''.join([area,'/assembly/','cp0',cp,'-/454NewblerMetrics.txt'])
        else:
            filename=''.join([area,'/assembly/','cp',cp,'-/454NewblerMetrics.txt'])
        with open(filename,'r') as In:
            number=re.
