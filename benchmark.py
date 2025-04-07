import time
import numpy as np
import random

d = 256
q = 3329  
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


# Function to compute a mod (X^{d/2} - r) and a mod (X^{d/2} + r)
def build_residue_tree(tree, node, mod_value, depth):
    if depth == 2:
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
    
def NTT(poly1, poly2, node, depth, list_):
    if depth == 2:
        new_poly = (poly1 * poly2) % ( X^depth + node.residue)
        list_.append(new_poly)
        return None
    else:
        n = depth // 2
        left = node.left
        right = node.right
        new_poly1 = (poly1) % ( X^depth + node.residue)
        new_poly2 = (poly2) % ( X^depth + node.residue)
        b = NTT(new_poly1, new_poly2, node.left, n, list_) 
        c = NTT(new_poly1, new_poly2, node.right, n, list_) 
        return None

def INTT(ntt_list, node, depth):
    if depth == 4:
        a = [0] * depth
        r = node.right
        r = R(r.residue)
        b, c = ntt_list
        for i in range(0, depth // 2):
            a[i] = (b[i] + c[i]) / 2
            a[i + depth // 2] = (r^(-1) * (b[i] - c[i])) / 2
        return PR(a)
    else:
        b = INTT(ntt_list[:depth // 4], node.left, depth // 2) 
        c = INTT(ntt_list[depth // 4:], node.right, depth // 2) 
        r = node.right
        r = r.residue
        a = [0] * depth  

        for i in range(0, depth // 2):
            a[i] = (b[i] + c[i]) / 2
            a[i + depth // 2] = (r^(-1) * (b[i] - c[i])) / 2

        return PR(a)  

def multiply(a, b, tree, dd):
    list_ntt = []
    NTT(a, b, tree.root, dd, list_ntt)
    intt = INTT(list_ntt, tree.root, dd)
    return intt

# Generar un polinomio f de grado d (256)
f = PR(list(np.random.randint(0, 2, d)))

# Medir tiempo para construir el árbol y realizar la multiplicación NTT
start_tree_build = time.time()
tree = ResidueTree(1)
build_residue_tree(tree, tree.root, q, depth=d)
end_tree_build = time.time()

start_multiply = time.time()
result = multiply(f, f, tree, d)
end_multiply = time.time()

# Medir el tiempo de la multiplicación normal
start_normal = time.time()
g = f * f
end_normal = time.time()

# Imprimir los resultados y los tiempos
print("Resultado NTT:", result)
print("Resultado normal:", Rk(g))
print("\nTiempo de construcción del árbol:", end_tree_build - start_tree_build, "segundos")
print("Tiempo de multiplicación NTT:", end_multiply - start_multiply, "segundos")
print("Tiempo de multiplicación normal:", end_normal - start_normal, "segundos")

# Función para la multiplicación manual
def multiplicacion_manual(p, q):
    resultado = 0
    for exp_p, coef_p in p.dict().items():
        for exp_q, coef_q in q.dict().items():
            resultado += coef_p * coef_q * x**(exp_p + exp_q)
    return resultado







d = 256
q = 3329  # Prime modulus
R = IntegerModRing(q)  # Anillo de enteros módulo q
PR.<X> = PolynomialRing(R)  # Anillo de polinomios en Z_q
Rk = PR.quotient(X^d + 1)  

def multiplicacion_manual(p, q):
    resultado = 0
    for exp_p, coef_p in p.dict().items():
        for exp_q, coef_q in q.dict().items():
            coef_p = R(coef_p)
            coef_q = R(coef_q)
            resultado += coef_p * coef_q * X**(exp_p + exp_q)
    return resultado


# Medir tiempo para la multiplicación manual
start_manual = time.time()
resultado_manual = Rk(multiplicacion_manual(f, f))
end_manual = time.time()

# Medir tiempo para la multiplicación normal
start_normal = time.time()
g = f * f
end_normal = time.time()

# Imprimir los resultados y los tiempos
print("Resultado manual:", resultado_manual)
print("Resultado normal:", Rk(g))
print("\nTiempo de multiplicación manual:", end_manual - start_manual, "segundos")
print("Tiempo de multiplicación normal:", end_normal - start_normal, "segundos")
