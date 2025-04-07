# SageMath implementation of polynomial multiplication using CRT in Z_q[X]/(X^d + α)
import numpy as np
import random
# Parámetros de Kyber-1024
q = 3329  # Módulo primo
R = IntegerModRing(q)  # Anillo de enteros módulo q
n = 256  # Dimensión del polinomio
PR.<X> = PolynomialRing(R)  # Anillo de polinomios en Z_q
Rk = PR.quotient(X^n + 1)  # Anillo cociente Z_q[X]/(X^n + 1)
k = 4  # Parámetro de Kyber
eta = 2  # Parámetro de ruido binomial
du, dv = 11, 5  # Parámetros de compresión

class ResidueTreeNode:
    def __init__(self, residue):
        self.residue = residue  # Residuo almacenado en el nodo
        self.left = None  # Hijo izquierdo
        self.right = None  # Hijo derecho

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



r_list=[]
def build_residue_tree(tree, node, mod_value, depth):
    if depth == 2:
        r_list.append(node.residue)
        return

    try:
        r1 = mod(-node.residue, mod_value).sqrt()  # r1 = sqrt(mod(-r, q))
        r_inv = R(r1)^-1  # r1^-1
        ri = (node.residue * R(r1)^-1) % mod_value  # ri = r * r1^-1

        tree.insert(node.residue, left_residue=ri, right_residue=r1)

        left_node = tree.find(tree.root, r1)
        right_node = tree.find(tree.root, ri)
        build_residue_tree(tree, left_node, mod_value, depth//2)
        build_residue_tree(tree, right_node, mod_value, depth//2)

    except ValueError:
        return
    
def NTT_prev(poly1, node, depth, list_):
    if depth == 2:
        new_poly = PR(poly1.list()) % ( X^depth + node.residue)
        list_.append(new_poly)
        return None
    else:
        d=depth//2
        left = node.left
        right = node.right
        new_poly1 = PR(poly1.list()) % ( X^depth + node.residue)
        b = NTT_prev(new_poly1, node.left, d, list_) 
        c = NTT_prev(new_poly1, node.right, d, list_) 
        return None

def INTT_prev(ntt_list, node, depth):
    if depth==4:
        a=[0]*depth
        r=node.right
        r=R(r.residue)
        b, c = ntt_list
        for i in range(0, depth//2):
            a[i] = (b[i] + c[i]) / 2
            a[i + depth//2] = (r^(-1) * (b[i] - c[i])) / 2
        return PR(a)
    else:
        b = INTT_prev(ntt_list[:depth//4], node.left, depth//2) 
        c = INTT_prev(ntt_list[depth//4:], node.right, depth//2) 
        r = node.right
        r=r.residue
        a = [0] * depth  

        for i in range(0, depth//2):
            a[i] = (b[i] + c[i]) / 2
            a[i + depth//2] = (r^(-1) * (b[i] - c[i])) / 2

        return PR(a)  
def NTT(poly):
    list_ = []
    NTT_prev(poly, tree.root, n, list_)
    return list_

def INTT(poly):
    result = INTT_prev(poly, tree.root, n)
    return result
def multiply(a, b):
    list_a = NTT(a)
    list_b = NTT(b)
    list_ntt=[]
    
    listt = r_list[::-1]
    
    for i in range(0, n//2):
        list_ntt.append((list_a[i]*list_b[i]) % (X^2+listt[i]))
    
    intt = INTT(list_ntt)
    return intt


def mod_switching(x, p, q1):
    x_ = [c.lift() for c in x.list()]
    f_ = p / q1
    m_bin_coeffs = [floor(c * f_ + 0.5) % p for c in x_]
    return Rk(m_bin_coeffs)

def binomial(eta):
    return sum(random.randint(0, 1) for _ in range(eta)) - sum(random.randint(0, 1) for _ in range(eta))

def binomial_sample(eta1):
    return Rk([binomial(eta1) for _ in range(n)])
tree = ResidueTree(1)

build_residue_tree(tree, tree.root, q, depth=n)

# Implementación de multiplicación manual usando bucles for
def multiply_vector_matrix(v, M):
    resultado = []
    for j in range(k): 
        suma = Rk(0)
        for i in range(k):  
            suma += multiply(v[i], M[j][i]) #v[i] * M[j, i]
        resultado.append(suma)

    return vector(Rk, resultado)

def multiply_vector_vector(v1, v2):
    p=sum([multiply(v1[i], v2[i]) for i in range(k)])
    return Rk(p)

# Key generation (sin cambios en la estructura)
def keygen():
    A = [[Rk([R(random.randint(0, q - 1)) for _ in range(n)]) for _ in range(k)] for _ in range(k)]
    s = [binomial_sample(eta) for _ in range(k)]
    e = [binomial_sample(eta) for _ in range(k)]
    
    #t = multiply_vector_matrix(s, A) + e
    s_A = multiply_vector_matrix(s, A)
    
    t = [s_A[i] + e[i] for i in range(k)] #multiply_vector_matrix(s_ntt, A_ntt) + e_ntt
    
    return (A, t), s

def encrypt(pk, m):
    A, t = pk
    r = [binomial_sample(eta) for _ in range(k)]
    e1 = [binomial_sample(eta) for _ in range(k)]
    e2 = binomial_sample(eta)
    A_t = [[A[j][i] for j in range(k)] for i in range(k)]
    r_A = multiply_vector_matrix(r, A_t)
    u = [r_A[i] + e1[i] for i in range(k)] #multiply_vector_matrix(r, A_t) + e1
    v =  multiply_vector_vector(r, t) + e2[0] + (q // 2) * Rk(m)
    
    u_compressed = [mod_switching(x, 2^du, q) for x in u]
    v_compressed = mod_switching(v, 2^dv, q)

    return u_compressed, v_compressed

def decrypt(sk, cipher):
    u, v = cipher  
    u_ = [mod_switching(x, q, 2^du) for x in u]
    v_ = mod_switching(v, q, 2^dv)

    m_dec = v_ - multiply_vector_vector(u_, sk)
    m_dec_coeffs = mod_switching(m_dec, 2, q)

    return Rk(m_dec_coeffs)

pk, sk = keygen()
m = list(np.random.randint(0, 2, n))
print(f"Mensaje original: {m}")

ciphertext = encrypt(pk, m)
m_decrypted = decrypt(sk, ciphertext)

print(f"Mensaje desencriptado: {m_decrypted}")

print("¿Coinciden el mensaje original y el desencriptado?:", m == [int(c) for c in m_decrypted.list()])