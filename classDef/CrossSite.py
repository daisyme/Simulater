#!/usr/bin/python
# -*- coding: UTF-8 -*-

class CrossSite(object):
    def __init__(self, parent, ancestor, location):
        self.pid = parent
        self.aid = ancestor
        self.loc = location
    def __repr__(self):
        return repr((self.pid, self.aid, self.loc))
    def __lt__(self, other):
        return self.loc < other.loc
    def __gt__(self, other):
        return self.loc > other.loc
    def __eq__(self, other):
        return self.loc == other.loc
    def __le__(self, other):
        return self.loc <= other.loc
    def __ge__(self, other):
        return self.loc >= other.loc
    def __ne__(self, other):
        return self.loc != other.loc
    def cut(self, other):
        return CrossSite(self.pid, other.aid, self.loc)
    def nextgen(self, pid):
        return CrossSite(pid, self.aid, self.loc)

#now the sorted is not within the structure which may slow down the speed
#TODO: improve the structure