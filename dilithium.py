#DILITHIUM
q = 2^23 - 2^13 + 1       
d = 256         
n, m = 6, 5
beta = 4       
gamma = 49*4     
boundary = beta + gamma  
delta_s = (q - 1)/32 - 1
R.<x> = PolynomialRing(GF(q))
f = x^d + 1   
Rq = R.quotient(f, 'a')

def key_gen():
    #TODO
    return 

def prover(A, t, s1, s2, c):
    #TODO
    return 

def verifier(A, t, w, c, response):
    #TODO
    return 

def simulate_protocol():
    #TODO
    return 

simulate_protocol()