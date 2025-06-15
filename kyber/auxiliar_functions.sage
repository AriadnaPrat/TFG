# Performs naive polynomial multiplication between p and q_
def trivial_multiplication(p, q_):
    result = 0
    for exp_p, coef_p in p.dict().items():
        for exp_q, coef_q in q_.dict().items():
            coef_p = R(coef_p)
            coef_q = R(coef_q)
            result += coef_p * coef_q * X**(exp_p + exp_q)
    return result

# Multiplies corresponding elements of two lists of polynomials and reduces each product modulo (X^2 + r_i)
def multiply(a, b, d):
    list_ntt=[]
    listt = r_list[::-1]
    for i in range(0, d//2):
        list_ntt.append(trivial_multiplication(a[i],b[i]) % (X^2+listt[i]))
    return list_ntt

# Switches modulus from q1 to p
def mod_switching(x, p, q1):
    x_ = [c.lift() for c in x.list()] # Get integer coefficients
    f_ = p / q1
    m_bin_coeffs = [floor(c * f_ + 0.5) % p for c in x_]
    return Rk(m_bin_coeffs)

# Samples from a centered binomial distribution with parameter eta
def binomial(eta):
    return sum(random.randint(0, 1) for _ in range(eta)) - sum(random.randint(0, 1) for _ in range(eta))

# Creates a polynomial with coefficients sampled from centered binomial distribution
def binomial_sample(eta, d):
    return PR([binomial(eta) for _ in range(d)])

# Multiplies a vector by a matrix (both of polynomials), returns resulting vector
def multiply_vector_matrix(v, M, dimension, d_):
    resultado = []
    for j in range(dimension): 
        suma = []
        for i in range(dimension):  
            suma.append(multiply(v[i], M[j][i], d_)) 
        sum_ = [sum([suma[j][i] for j in range(dimension)]) for i in range(d_//2)]
        resultado.append(sum_)
    return resultado

# Multiplies two vectors of polynomials
def multiply_vector_vector(v1, v2, dimension, d_):
    result = [multiply(v1[i], v2[i], n) for i in range(dimension)]
    p=[sum([result[i][j] for i in range(dimension)]) for j in range(d_//2)]
    return p

# Generates a list of random degree-1 polynomials of the form a + bX
def generate_ntt(d_):
    return [R(random.randint(0, q - 1)) + R(random.randint(0, q - 1)) * X for _ in range(d_//2)]

# Generates a matrix of random polynomials to be used with NTT operations
def create_matrix_ntt(d_, dimension):
    return [[generate_ntt(d_) for _ in range(dimension)] for _ in range(dimension)]

# Creates a vector of polynomials with coefficients from centered binomial distribution
def create_vector(eta, dimension):
    return [binomial_sample(eta, n) for _ in range(dimension)]

def sum_vectors(a, b, dimension):
    return [a[i] + b[i] for i in range(dimension)] 

# Transposes a square matrix
def matrix_traspose(matrix, dimension):
    return [[matrix[j][i] for j in range(dimension)] for i in range(dimension)]

def NTT_vector(vector, d):
    return [NTT(i, d) for i in vector]

def INTT_vector(vector):
    return [INTT(i) for i in vector]

# Converts a byte array to a polynomial with coefficients from centered binomial distribution
def cbd(data, eta):
    bits = ''.join(f'{b:08b}' for b in data)
    poly = []
    for i in range(0, n): 
        a = sum(int(bits[2*i*eta + j]) for j in range(eta))
        b = sum(int(bits[2*i*eta + eta + j]) for j in range(eta))
        poly.append(a - b)
    return PR(poly)

# Pseudorandom function based on SHAKE-256; returns byte string of length 64*eta
def PRF(seed, nonce, eta):
    total = seed + format(nonce, 'b')
    x = hashlib.shake_256(total.encode()).digest(64*eta)
    return x

# This function Hash processes message m and public key pk
def H(m, pk):
    A, t = pk 
    A_ = sum(sum(A, []), [])
    t_ = sum(t, [])
    total = A_ + t_
    coef_binary = [[bin(j) for j in i.list()] for i in total]
    coef_binary = sum(coef_binary, [])
    m_binary = ''.join(coef_binary)
    mm = ''.join([bin(j) for j in m])
    hash_obj = hashlib.sha3_256(m_binary.encode()).digest()
    tt = ''.join([mm, ''.join(format(byte, '08b') for byte in hash_obj)]) 
    hash_ = hashlib.sha3_512(tt.encode()).digest()
    hash_ = ''.join(f'{byte:08b}' for byte in hash_)
    return hash_

# Generates a vector of polynomials deterministically from a seed using PRF and cbd
def create_vector_(eta , d , seed ):
    vec = []
    nonce = 0
    for _ in range(d):
        output = PRF(seed, nonce, eta)
        poly = cbd(output, eta)
        vec.append(poly)
        nonce += 1
    return vec
