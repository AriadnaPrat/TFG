def trivial_multiplication(p, q_):
    result = 0
    for exp_p, coef_p in p.dict().items():
        for exp_q, coef_q in q_.dict().items():
            coef_p = R(coef_p)
            coef_q = R(coef_q)
            result += coef_p * coef_q * X**(exp_p + exp_q)
    return result
    
def multiply(a, b, d):
    list_ntt=[]
    listt = r_list[::-1]
    for i in range(0, d//2):
        list_ntt.append(trivial_multiplication(a[i],b[i]) % (X^2+listt[i]))
    return list_ntt


def mod_switching(x, p, q1):
    x_ = [c.lift() for c in x.list()]
    f_ = p / q1
    m_bin_coeffs = [floor(c * f_ + 0.5) % p for c in x_]
    return Rk(m_bin_coeffs)

def binomial(eta):
    return sum(random.randint(0, 1) for _ in range(eta)) - sum(random.randint(0, 1) for _ in range(eta))

def binomial_sample(eta, d):
    return Rk([binomial(eta) for _ in range(d)])

def multiply_vector_matrix(v, M, dimension, d_):
    resultado = []
    for j in range(dimension): 
        suma = []
        for i in range(dimension):  
            suma.append(multiply(v[i], M[j][i], d_)) 
        sum_ = [sum([suma[j][i] for j in range(dimension)]) for i in range(d_//2)]
        resultado.append(sum_)
    return resultado

def multiply_vector_vector(v1, v2, dimension, d_):
    result = [multiply(v1[i], v2[i], n) for i in range(dimension)]
    p=[sum([result[i][j] for i in range(dimension)]) for j in range(d_//2)]
    return p

def generate_ntt(d_):
    return [R(random.randint(0, q - 1)) + R(random.randint(0, q - 1)) * X for _ in range(d_//2)]

def create_matrix_ntt(d_, dimension):
    return [[generate_ntt(d_) for _ in range(dimension)] for _ in range(dimension)]

def create_vector(eta, dimension):
    return [binomial_sample(eta, n) for _ in range(dimension)]

def sum_vectors(a, b, dimension):
    return [a[i] + b[i] for i in range(dimension)] 

def matrix_traspose(matrix, dimension):
    return [[matrix[j][i] for j in range(dimension)] for i in range(dimension)]

def NTT_vector(vector, d):
    return [NTT(i, d) for i in vector]

def INTT_vector(vector):
    return [INTT(i) for i in vector]


def cbd(data, eta):
    bits = ''.join(f'{b:08b}' for b in data)
    poly = []
    for i in range(0, n): 
        a = sum(int(bits[2*i*eta + j]) for j in range(eta))
        b = sum(int(bits[2*i*eta + eta + j]) for j in range(eta))
        poly.append(a - b)
    return PR(poly)

def PRF(seed, nonce, eta):
    total = seed + format(nonce, 'b')
    x = hashlib.shake_256(total.encode()).digest(64*eta)
    return x

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

def create_vector_(eta , d , seed ):
    vec = []
    nonce = 0
    for _ in range(d):
        output = PRF(seed, nonce, eta)
        poly = cbd(output, eta)
        vec.append(poly)
        nonce += 1
    return vec