#DILITHIUM
import random
import numpy as np
import time
q =  2^23 - 2^13 + 1 
d = 256 #256         
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
    list_a = NTT(a)
    list_b = NTT(b)
    list_ntt=[]
    
    listt = r_list[::-1]
    for i in range(d):
        list_ntt.append(list_a[i]*list_b[i])
    
    intt = INTT(list_ntt)
    return intt

f = PR(list(np.random.randint(0, 2, d)))

start_tree_build = time.time()
tree, r_list = build_tree(1, d)
end_tree_build = time.time()

start_multiply = time.time()
result = multiply(f, f)
end_multiply = time.time()

start_normal = time.time()
sage_f = Rk(f * f)
end_normal = time.time()


def trivial_multiplication(p, q):
    result = 0
    for exp_p, coef_p in p.dict().items():
        for exp_q, coef_q in q.dict().items():
            coef_p = R(coef_p)
            coef_q = R(coef_q)
            result += coef_p * coef_q * X**(exp_p + exp_q)
    return result


start_manual = time.time()
resultado_manual = Rk(trivial_multiplication(f, f))
end_manual = time.time()
#print("NTT result:", result)
#print("SageMath Result:", sage_f)
#print("\nTrivial Multiplication:", resultado_manual)
print("¿The multiplications got the same result?:", resultado_manual == result == sage_f)
print("Time construction Tree:", end_tree_build - start_tree_build, "segundos")
print("NTT Complexity Time:", end_multiply - start_multiply, "segundos")
print("SageMath Algebra Complexity Time:", end_normal - start_normal, "segundos")
print("Trivial Multiplication:", end_manual - start_manual, "segundos")
