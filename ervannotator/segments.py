#######################################################################
## Author: 	Thanos Mourikis
## Date: 	2014
## Purpose:	Segments class for overlap elimination
#######################################################################


from itertools import izip
from bisect import bisect_left


class Segments: ## Class to resolve overlapping regions

    def __init__(self, min_hsp_distance):
        self.Ls = []
        self.Rs = []
        self.min_hsp_distance = min_hsp_distance
        
    def add(self, l, r):
        ## Instead of throwing an error reverse the coords and pass them to object
        if (l <= r):            
            i = bisect_left(self.Ls, l)
            self.Ls.insert(i, l)
            self.Rs.insert(i, r)
        else:
            i = bisect_left(self.Ls, l)
            self.Ls.insert(i, r)
            self.Rs.insert(i, l)

    def __iter__(self):
        if len(self.Ls) == 0:
            raise StopIteration
        l = self.Ls[0]
        r = self.Rs[0]
        for _l, _r in izip(self.Ls[1:], self.Rs[1:]):
            if _l <= r + self.min_hsp_distance: ## This line can be altered according to buffer
                r = max(r, _r)
            else:
                yield (l, r)
                l = _l
                r = _r
        yield (l, r)
        
    def compact(self):
        _Ls = []
        _Rs = []
        for (l, r) in self:
            _Ls.append(l)
            _Rs.append(r)
        self.Ls = _Ls
        self.Rs = _Rs
        
    def wipe(self):
        del self.Ls[:]
        del self.Rs[:]
    
    def __repr__(self):
        return str((self.Ls, self.Rs))
    
    def __str__(self):
        return str(list(self))

    def __len__(self):
        self.compact()
        n = len(self.Ls)
        assert(n == len(self.Rs))
        return n