from legsob import *
import sys

if len(sys.argv) != 3:
    print "Usage: legsobpolys.py mu order"
    exit()
    
mu, n = tuple(sys.argv[1:])
mu = float(mu)
n = int(n)



lsp = LegendreSobolevPolynomials()
lsp.mu = [1,mu]
lsp.n = n
print lsp.n
print lsp.mu[1]
for p in lsp.LSPdiffN:
    s = ""
    for c in p:
        s += str(round(c, 17)) + " "
    print s
s = ""
for c in lsp.e0N:
    s += str(round(c, 17)) + " "
print s
s = ""
for c in lsp.e1N:
    s += str(round(c, 17)) + " "
print s
for p in lsp.LSPN:
    s = ""
    for c in p:
        s += str(round(c, 17)) + " "
    print s
