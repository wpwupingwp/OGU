from multiprocessing.dummy import Pool as thread_pool
from multiprocessing import call

urls=[
        'ftp://Novo_NH150704_100909:MGYzNWEx@ftpdata.novogene.cn:2300/results/RawData/C16_1.fq.gz',
        'ftp://Novo_NH150704_100909:MGYzNWEx@ftpdata.novogene.cn:2300/results/RawData/C16_2.fq.gz',
        'ftp://Novo_NH150704_100909:MGYzNWEx@ftpdata.novogene.cn:2300/results/RawData/C17_1.fq.gz',
        'ftp://Novo_NH150704_100909:MGYzNWEx@ftpdata.novogene.cn:2300/results/RawData/C17_2.fq.gz',
        'ftp://Novo_NH150704_100909:MGYzNWEx@ftpdata.novogene.cn:2300/results/RawData/C18_1.fq.gz',
        'ftp://Novo_NH150704_100909:MGYzNWEx@ftpdata.novogene.cn:2300/results/RawData/C18_2.fq.gz',
        'ftp://Novo_NH150704_100909:MGYzNWEx@ftpdata.novogene.cn:2300/results/RawData/C19_1.fq.gz',
        'ftp://Novo_NH150704_100909:MGYzNWEx@ftpdata.novogene.cn:2300/results/RawData/C19_2.fq.gz'
        ]

pool = thread_pool(8)
results = pool.map(down, urls)

def down(url):
    return call(['wget', url])
