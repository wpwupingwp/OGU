from Bio.Blast.Applications import NcbiblastnCommandline as nb
from multiprocessing import cpu_count


def blast(query_file, db_file, result='BLASTResult.xml'):
    """Here we use "max_hsps" to restrict only the first hsp, use
    "max_target_seqs" to restrict only the first matched sequence.
    """
    cmd = nb(num_threads=cpu_count(),
             query=query_file,
             db=db_file,
             task='blastn',
             evalue=arg.evalue,
             max_hsps=1,
             max_target_seqs=1,
             outfmt=5,
             out=result)
    stdout, stderr = cmd()
    return result
