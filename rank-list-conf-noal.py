import sys
import os
from time import *

if len(sys.argv) != 3:
    print "Usage: rank-list-conf-noal.py file-with-test-samples file-with-rank-lists"
    print "Determines the rank of the correct class for each test sample."
    exit()

labelfn = sys.argv[1]
r1fn = sys.argv[2]

flabels_file = open(labelfn, "rt")
flabels = flabels_file.readlines()
cls = {}
rk = {}
classname = {}
relevant = set([])
for ln in flabels:
    ch, cl = tuple(ln[:-1].split(" "))
    if cl not in cls:
        cls[cl] = []    
    cls[cl].append(ch)
    rk[ch] = []
    classname[ch] = cl
    relevant.add(cl)
    
flabels_file.close()

ln = open(r1fn,"rt").readlines()

for i in range(len(ln)):
    u = ln[i].split()
    ch = u[0]
    orig = set(classname[ch].split("-")[:-1])

    rl = u[1:]
    mink = len(rl) / 2
    for k in range(len(rl) / 2):
        cur = set(rl[2*k].split("-")[:-1])
        if orig.intersection(cur):
            mink = min(mink, k)
            break
    print ch, mink+1, rl[1]

