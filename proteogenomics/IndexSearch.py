import _IndexSearch as _c_impl

class IndexSearch(object):
    """Search for protein string using Seqan k-mer index"""

    def __init__(self,dbStr):
        self._c_obj = _c_impl.new_IndexSearch(dbStr)

    def find(self,pepStr):
        return _c_impl.IndexSearch_find(self._c_obj,pepStr)


if __name__ == "__main__":
    dbStr = "".join("""\
SATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKN
DMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGK
VQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGP
CTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIR
IQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNC
GGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSN
ITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKR""".\
    strip().split())
    pepStr = "VQKEYAFFY"

    for i in xrange(10**6):
        ins_res = IndexSearch(dbStr).find(pepStr)
        find_res = [dbStr.find(pepStr)]
        #print "IndexSearch results:  ", ins_res
        #print "Python Search results:", find_res
        assert ins_res == find_res


    
