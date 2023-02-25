# collect third-party files from original url
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen

urls = [
    'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-win64.tar.gz',
    'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-macosx.tar.gz',
    'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz',
    'https://github.com/iqtree/iqtree2/releases/download/v2.2.2.2/iqtree-2.2.2.2-Linux.tar.gz',
    'https://github.com/iqtree/iqtree2/releases/download/v2.2.2.2/iqtree-2.2.2.2-MacOSX.zip',
    'https://github.com/iqtree/iqtree2/releases/download/v2.2.2.2/iqtree-2.2.2.2-Windows.zip',
    'https://mafft.cbrc.jp/alignment/software/mafft-7.511-win64-signed.zip',
    'https://mafft.cbrc.jp/alignment/software/mafft-7.511-mac.zip?signed',
    'https://mafft.cbrc.jp/alignment/software/mafft-7.511-linux.tgz',
]


def get(url) -> str:
    # return url for fail
    filename = url.split('/')[-1].replace('?signed', '')
    print('Fetching', url)
    try:
        down = urlopen(url, timeout=100)
    except Exception:
        return url
    with open(filename, 'wb') as out:
        out.write(down.read())
    return filename


with ThreadPoolExecutor() as pool:
    for url in urls:
        pool.submit(get, url)