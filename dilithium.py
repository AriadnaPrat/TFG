#DILITHIUM
q = 2^23 - 2^13 + 1       
d = 256         
n, m = 6, 5
beta = 4       
gamma = 49*4     
beta_prima = 2^19 - gamma + 1
delta_s = (q - 1)/32 - 1
R.<x> = PolynomialRing(GF(q))
f = x^d + 1   
Rq = R.quotient(f, 'a')

def binomial(eta):
    return sum(random.randint(0, 1) for _ in range(eta)) - sum(random.randint(0, 1) for _ in range(eta))

def binomial_sample(eta):
    return Rk([binomial(eta) for _ in range(n)])

def key_gen():
    #TODO
    return 
def create_vector(eta, d):
    return [binomial_sample(eta) for _ in range(d)]

def multiply_vector_matrix(v, M):
    resultado = []
    for j in range(k): 
        suma = []
        for i in range(k):  
            suma.append(multiply(v[i], M[j][i])) #v[i] * M[j, i]
        sum_ = [sum([suma[j][i] for j in range(k)]) for i in range(n//2)]
        resultado.append(sum_)
    return resultado

def prover(A, t, s1, s2, c):
    y_1 = create_vector(gamma + beta_prima, m)
    y_2 = create_vector(gamma + beta_prima, n)
    
    y_1_A = multiply_vector_matrix(y_1, A)
    
    
    #TODO
    return 

def verifier(A, t, w, c, response):
    #TODO
    return 

def simulate_protocol():
    #TODO
    return 

simulate_protocol()