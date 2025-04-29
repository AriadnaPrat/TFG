# SageMath implementation of polynomial multiplication using CRT in Z_q[X]/(X^d + α)
import numpy as np
import random
import time

# Parámetros de Kyber-1024
q = 3329  
R = IntegerModRing(q)  
n = 256  
PR.<X> = PolynomialRing(R)  
Rk = PR.quotient(X^n + 1)  
k = 4  
eta1, eta2 = 2, 2
du, dv = 11, 5  

class ResidueTreeNode:
    def __init__(self, residue):
        self.residue = residue  
        self.left = None  
        self.right = None  

class ResidueTree:
    def __init__(self, root_residue):
        self.root = ResidueTreeNode(root_residue)

    def insert(self, parent_residue, left_residue=None, right_residue=None):
        parent_node = self.find(self.root, parent_residue)
        if parent_node:
            if left_residue is not None:
                parent_node.left = ResidueTreeNode(left_residue)
            if right_residue is not None:
                parent_node.right = ResidueTreeNode(right_residue)

    def find(self, node, residue):
        if node is None:
            return None
        if node.residue == residue:
            return node
        return self.find(node.left, residue) or self.find(node.right, residue)

    def display(self, node=None, level=0):
        if node is None:
            node = self.root
        print("  " * level + f"({node.residue})")
        if node.left:
            self.display(node.left, level + 1)
        if node.right:
            self.display(node.right, level + 1)

def build_residue_tree(tree, node, mod_value, depth, r_list):
    if depth == 2:
        r_list.append(node.residue)
        return

    try:
        r1 = mod(-node.residue, mod_value).sqrt()  
        r_inv = R(r1)^-1  # r1^-1
        ri = (node.residue * r_inv) % mod_value  
        tree.insert(node.residue, left_residue=ri, right_residue=r1)
        left_node = tree.find(tree.root, r1)
        right_node = tree.find(tree.root, ri)
        build_residue_tree(tree, left_node, mod_value, depth//2, r_list)
        build_residue_tree(tree, right_node, mod_value, depth//2, r_list)

    except ValueError:
        return
    
def build_tree(alpha, degree):
    r_list = []
    tree = ResidueTree(alpha)
    build_residue_tree(tree, tree.root, q, degree, r_list)
    return tree, r_list

def NTT_prev(a, node, depth, list_):
    if depth == 2:
        list_.append(a[0]+a[1]*X)
        return None
    else:
        d = depth // 2
        r=node.right
        r=R(r.residue)
        b = [0]*d
        c = [0]*d
        poly = a
        for i in range(d):
            b[i] = (poly[i] + r *poly[i + d])
            c[i] = (poly[i] - r *poly[i + d])
            
        NTT_prev(b, node.left, d, list_) 
        NTT_prev(c, node.right, d, list_) 
        return None

def INTT_prev(ntt_list, node, depth):
    a=[0]*depth
    r=node.right
    r=R(r.residue)
    d = depth//2
    
    if depth==4:
        b, c = ntt_list
        degree_half = n // 2
        for i in range(0, d):
            a[i] = (b[i] + c[i]) / degree_half
            a[i + depth//2] = (r^(-1) * (b[i] - c[i])) / degree_half
        return PR(a)
    else:
        b = INTT_prev(ntt_list[:d//2], node.left, depth//2) 
        c = INTT_prev(ntt_list[d//2:], node.right, depth//2) 
        for i in range(0, d):
            a[i] = (b[i] + c[i])
            a[i + d] = (r^(-1) * (b[i] - c[i]))

        return PR(a)  
    
def NTT(poly):
    list_ = []
    NTT_prev(poly, tree.root, n, list_)
    return list_

def INTT(poly):
    result = INTT_prev(poly, tree.root, n)
    return result

def trivial_multiplication(p, q):
    resultado = 0
    for exp_p, coef_p in p.dict().items():
        for exp_q, coef_q in q.dict().items():
            coef_p = R(coef_p)
            coef_q = R(coef_q)
            resultado += coef_p * coef_q * X**(exp_p + exp_q)
    return resultado
    
def multiply(a, b):
    list_ntt=[]
    listt = r_list[::-1]
    for i in range(0, n//2):
        list_ntt.append(trivial_multiplication(a[i],b[i]) % (X^2+listt[i]))
    
    return list_ntt


def mod_switching(x, p, q1):
    x_ = [c.lift() for c in x.list()]
    f_ = p / q1
    m_bin_coeffs = [floor(c * f_ + 0.5) % p for c in x_]
    return Rk(m_bin_coeffs)

def binomial(eta):
    return sum(random.randint(0, 1) for _ in range(eta)) - sum(random.randint(0, 1) for _ in range(eta))

def binomial_sample(eta):
    return Rk([binomial(eta) for _ in range(n)])

def multiply_vector_matrix(v, M):
    resultado = []
    for j in range(k): 
        suma = []
        for i in range(k):  
            suma.append(multiply(v[i], M[j][i])) 
        sum_ = [sum([suma[j][i] for j in range(k)]) for i in range(n//2)]
        resultado.append(sum_)
    return resultado

def multiply_vector_vector(v1, v2):
    result = [multiply(v1[i], v2[i]) for i in range(k)]
    p=[sum([result[i][j] for i in range(k)]) for j in range(n//2)]
    return p

def generate_ntt():
    return [R(random.randint(0, q - 1)) + R(random.randint(0, q - 1)) * X for _ in range(n//2)]

def create_matrix_ntt():
    return [[generate_ntt() for _ in range(k)] for _ in range(k)]

def create_vector(eta, dimension):
    return [binomial_sample(eta) for _ in range(dimension)]

def sum_vectors(a, b):
    return [a[i] + b[i] for i in range(k)] 

def matrix_traspose(matrix):
    return [[matrix[j][i] for j in range(k)] for i in range(k)]

def NTT_vector(vector):
    return [NTT(i) for i in vector]

def INTT_vector(vector):
    return [INTT(i) for i in vector]



def keygen():
    A_ntt = create_matrix_ntt()
    s = create_vector(eta1, k)
    e = create_vector(eta1, k) 
    s_ntt = NTT_vector(s)
    e_ntt = NTT_vector(e)
    s_A_ntt = multiply_vector_matrix(s_ntt, A_ntt)
    
    t_ntt = sum_vectors(s_A_ntt, e_ntt) 
    
    return (A_ntt, t_ntt), s_ntt


def cbd(data: bytes, eta: int) -> list:
    bits = ''.join(f'{b:08b}' for b in data)
    poly = []
    for i in range(0, n): #while len(poly) < n and i + 2 * eta <= len(data):
        a = sum(int(bits[2*i*eta + j]) for j in range(eta))
        b = sum(int(bits[2*i*eta + eta + j]) for j in range(eta))
        poly.append(a - b)
    return PR(poly)

def PRF(seed, nonce, eta):
    total = seed + format(nonce, 'b')
    x = hashlib.shake_256(total.encode()).digest(64*eta)
    return x
    
def create_vector_(eta: int, d: int, seed: bytes) -> tuple:
    vec = []
    nonce = 0
    for _ in range(d):
        output = PRF(seed, nonce, eta)
        poly = cbd(output, eta)
        vec.append(poly)
        nonce += 1
    return vec


def encrypt(pk, m, rrho = 0):
    A, t = pk
    r = create_vector_(eta1, k, rrho)
    e1 = create_vector_(eta2, k, rrho) 
    e2 = create_vector_(eta2, 1, rrho)[0] 
    r_ntt = NTT_vector(r) 
    A_t = matrix_traspose(A)
    r_A = multiply_vector_matrix(r_ntt, A_t)
    r_A_intt = INTT_vector(r_A) 
    u = sum_vectors(r_A_intt, e1) 
    r_t_intt = INTT(multiply_vector_vector(r_ntt, t))
    v =  r_t_intt + e2 + (q // 2) * Rk(m)
    
    u_compressed = [mod_switching(x, 2^du, q) for x in u]
    v_compressed = mod_switching(v, 2^dv, q)

    return u_compressed, v_compressed

def decrypt(sk, cipher):
    u, v = cipher  
    u_ = [mod_switching(x, q, 2^du) for x in u]
    v_ = mod_switching(v, q, 2^dv)
    u_ntt = NTT_vector(u_) 
    u_sk = multiply_vector_vector(u_ntt, sk)
    u_sk_intt = INTT(u_sk)
    m_dec = v_ - u_sk_intt
    m_dec_coeffs = mod_switching(m_dec, 2, q)

    return Rk(m_dec_coeffs)

import hashlib

def H(m, pk):
    A, t = pk 
    A_ = sum(sum(A, []), [])
    t_ = sum(t, [])
    total = A_ + t_
    coef_binarios = [[bin(j) for j in i.list()] for i in total]
    coef_binarios = sum(coef_binarios, [])
    mensaje_binario = ''.join(coef_binarios)
    mm = ''.join([bin(j) for j in m])
    hash_obj = hashlib.sha3_256(mensaje_binario.encode()).digest()
    tt = ''.join([mm, ''.join(format(byte, '08b') for byte in hash_obj)]) 
    hash_ = hashlib.sha3_512(tt.encode()).digest()
    hash_ = ''.join(f'{byte:08b}' for byte in hash_)
    return hash_
      
def KEM_keygen():
    pk, sk = keygen()
    return pk, sk

def KEM_encaps(pk):
    m = list(np.random.randint(0, 2, n))
    hash_ = H(m, pk)
    K, rrho = hash_[0:256], hash_[257:512]
    c = encrypt(pk, m, rrho)
    return K, c

def KEM_decaps(sk, c, K):
    m_prima = decrypt(sk, c)
    hash_ = H(m_prima, pk)
    K_prima, rrho_prima = hash_[0:256], hash_[257:512]
    c_prima = encrypt(pk, m_prima, rrho_prima)
    if c != c_prima:
        return None
    else:
        return K_prima

t0 = time.time()
tree, r_list = build_tree(1, n)
#pk, sk = keygen()
#m = list(np.random.randint(0, 2, n))

#ciphertext = encrypt(pk, m)
#m_decrypted = decrypt(sk, ciphertext)
#tiempo_total = time.time() - t0
#print(f"Original message: {m}")
#print(f"Decrypted message: {m_decrypted}")
#print("Does the original message match the decrypted one?:", m == [int(c) for c in m_decrypted.list()])
#print(f"Elapsed real time: {tiempo_total:.6f} seconds")

pk, sk = KEM_keygen()
K, c = KEM_encaps(pk)
k_shared=KEM_decaps(sk, c, K)
print("Keys:", pk, sk)
print("Encrypted:", c)
print("k_shared:", k_shared)
tiempo_total = time.time() - t0
print(f"Elapsed real time: {tiempo_total:.6f} seconds")