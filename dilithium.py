#DILITHIUM
import random
from random import choice, shuffle
import numpy as np
q =  2^23 - 2^13 + 1     
k = 4
d = 256 #256         
n, m = 6, 5
beta = 4       
gamma = 49*4     
beta_prima = 2^19 - gamma + 1
delta_s = (q - 1)/32 - 1
R = IntegerModRing(q)  
PR.<X> = PolynomialRing(R)  
Rk = PR.quotient(X^d + 1)  


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
    if depth == 1:
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

tree, r_list = build_tree(1, d)
def NTT_prev(a, node, depth, list_):
    if depth == 1:
        list_.append(a[0])
        return None
    else:
        dd = depth // 2
        r=node.right
        r=R(r.residue)
        b = [0]*dd
        c = [0]*dd
        for i in range(dd):
            b[i] = (a[i] + r *a[i + dd])
            c[i] = (a[i] - r *a[i + dd])
            
        NTT_prev(b, node.left, dd, list_) 
        NTT_prev(c, node.right, dd, list_) 
        return None

def INTT_prev(ntt_list, node, depth):
    a=[0]*depth
    dd = depth//2
    r=node.right
    r=R(r.residue)
    if depth==2:
        b, c = ntt_list
        a[0] = (b + c) / d
        a[1] = (r^(-1) * (b - c)) / d
        return PR(a)
    else:
        
        b = INTT_prev(ntt_list[:dd], node.left, dd) 
        c = INTT_prev(ntt_list[dd:], node.right, dd) 
        for i in range(dd):
            a[i] = (b[i] + c[i])
            a[i + dd] = (r^(-1) * (b[i] - c[i]))

        return PR(a)  
    
def NTT(poly):
    list_ = []
    NTT_prev(poly, tree.root, d, list_)
    return list_

def INTT(poly):
    result = INTT_prev(poly, tree.root, d)
    return result
    
def multiply(a, b):
    list_ntt=[]
    for i in range(d):
        list_ntt.append(a[i]*b[i])
    return list_ntt


def binomial(eta):
    return random.randint(-eta, eta) 

def binomial_sample(eta):
    return Rk([binomial(eta) for _ in range(d)])


def create_vector(eta, dimension):
    return [binomial_sample(eta) for _ in range(dimension)]

def sum_vectors(a, b, d1):
    return [a[i] + b[i] for i in range(d1)] 

def multiply_vector_matrix(v, M, d1, d2):
    resultado = []
    for j in range(d1): 
        suma = []
        for i in range(d2):  
            suma.append(multiply(v[i], M[j][i])) #v[i] * M[j, i]
        sum_ = [sum([suma[j][i] for j in range(m)]) for i in range(d)]
        resultado.append(sum_)
    return resultado

def generate_c():
    coef = [choice([-1, 1]) for _ in range(49)] + [0] * (d-49)

    shuffle(coef)
    
    return coef

def multiply_constant_vector(c, v):
    return [multiply(c,i) for i in v]

def create_matrix(d1, d2):
    return [[Rk([R(random.randint(0, q - 1)) for _ in range(d)]) for _ in range(d2)] for _ in range(d1)]

def NTT_vector(vector):
    return [NTT(i) for i in vector]

def INTT_vector(vector):
    return [INTT(i) for i in vector]

def is_in_beta_range(v, beta_bar):
    return all(-beta_bar <= x <= beta_bar for x in v)

class Prover:
    def __init__(self):
        self.y_1 = None
        self.y_2 = None
        self.s1 = None
        self.s2 = None
        self.z1 = None
        self.z2 = None
        self.A = None
        self.t = None
        
    def keygen(self):
        self.s1 = create_vector(beta, m)
        self.s2 = create_vector(beta, n)
        self.A = create_matrix(n, m)
        A_ntt = [NTT_vector(self.A[i]) for i in range(n)]
        s1_ntt = NTT_vector(self.s1)
        t_ = multiply_vector_matrix(s1_ntt, A_ntt, n, m)
        self.t = sum_vectors(INTT_vector(t_), self.s2, n)

        return self.A, self.t
        
    def step1(self):
        self.y_1 = create_vector(gamma + beta_prima, m)
        self.y_2 = create_vector(gamma + beta_prima, n)
        print(self.y_1)
        print(self.y_2)
        y_1_ntt = NTT_vector(self.y_1)
        A_ntt = [NTT_vector(self.A[i]) for i in range(n)]
        y_1_A_ntt = multiply_vector_matrix(y_1_ntt, A_ntt, n, m)
        y_1_A_intt = INTT_vector(y_1_A_ntt)
        self.w = sum_vectors(y_1_A_intt, self.y_2, n)
        return self.w
    
    def step2(self, c):
        s1_ = NTT_vector(self.s1)
        s2_ = NTT_vector(self.s2)
        c_ntt = NTT(c)
        z_1_ = multiply_constant_vector(c_ntt, s1_)
        z_2_ = multiply_constant_vector(c_ntt, s2_)
        
        z_1 = sum_vectors(INTT_vector(z_1_),  self.y_1, m)
        z_2 = sum_vectors(INTT_vector(z_2_), self.y_2, n)
        print("ddd:", z_1, z_2)
        if not is_in_beta_range(z_1, beta_prima) or not is_in_beta_range(z_2, beta_prima):
            return None, None
        else:
            return z_1, z_2
        
class Verifier:
    def __init__(self):
        self.c = None
        
    def step1(self, A, t):
        ones_neg_ones = [choice([-1, 1]) for _ in range(49)]
        zeros = [0] * (d - 49)
        self.c = ones_neg_ones + zeros
        shuffle(self.c)
        return PR(self.c)
    
    def step2():
        #TODO
        return


def simulate_protocol():
    v = Verifier()
    p = Prover()
    A, t = p.keygen()
    w = p.step1()
    c = v.step1(A, t)
    z1, z2 = p.step2(c)
    print(c)
    print(w)
    print(z1, z2)
    #TODO
    return 

simulate_protocol()