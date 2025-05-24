import numpy as np
import random
import time
import hashlib
import sys

# Kyber-1024 parameters

q, n = 3329, 256
R = IntegerModRing(q)  
PR.<X> = PolynomialRing(R)  
Rk = PR.quotient(X^n + 1)  
load("class_tree.sage") 
load("NTT_functions.sage")
load("auxiliar_functions.sage")
load("kyber.sage")
tree, r_list = build_tree(1, n, q)

def size_publickey(k_):
    bits_per_coef = ceil(log(q, 2))
    num_polynomials = k_^2 + k_
    total_coefficients = num_polynomials * n
    total_bits = total_coefficients * bits_per_coef
    total_bytes = ceil(total_bits / 8)
    return total_bits, total_bytes

def size_ct(k_):
    bits_per_coef = ceil(log(q, 2))
    num_polynomials = k_ + 1
    total_coefficients = num_polynomials * n
    total_bits = total_coefficients * bits_per_coef
    total_bytes = ceil(total_bits / 8)
    return total_bits, total_bytes

############################ kyber-512 ############################
t0 = time.time()
k = 2
kyber = Kyber_KEM(k, n, 3, 2, 10, 4, q)
pk, sk = kyber.KEM_keygen()
K, c = kyber.KEM_encaps(pk)
k_shared=kyber.KEM_decaps(sk, c, K)
t_total = time.time() - t0
print("k_shared:", k_shared)
print(f"Elapsed real time: {t_total:.6f} seconds")

total_bits_p, total_bytes_p = size_publickey(k)
total_bits_c, total_bytes_c = size_ct(k)
print(f"Public Key Size: {total_bits_p} bits ({total_bytes_p} bytes)")
print(f"Ciphertext: {total_bits_c} bits ({total_bytes_c} bytes)")

############################ kyber-768 ############################
t0 = time.time()
k = 3
kyber = Kyber_KEM(k, n, 2, 2, 10, 4, q)
pk, sk = kyber.KEM_keygen()
K, c = kyber.KEM_encaps(pk)
k_shared=kyber.KEM_decaps(sk, c, K)
print("k_shared:", k_shared)
t_total = time.time() - t0

print(f"Elapsed real time: {t_total:.6f} seconds")

total_bits_p, total_bytes_p = size_publickey(k)
total_bits_c, total_bytes_c = size_ct(k)
print(f"Public Key Size: {total_bits_p} bits ({total_bytes_p} bytes)")
print(f"Ciphertext: {total_bits_c} bits ({total_bytes_c} bytes)")

############################ kyber-1024 ############################
t0 = time.time()
k = 4
kyber = Kyber_KEM(k, n, 2, 2, 11, 5, q)
pk, sk = kyber.KEM_keygen()
K, c = kyber.KEM_encaps(pk)
k_shared=kyber.KEM_decaps(sk, c, K)
print("k_shared:", k_shared)
t_total = time.time() - t0

print(f"Elapsed real time: {t_total:.6f} seconds")

total_bits_p, total_bytes_p = size_publickey(k)
total_bits_c, total_bytes_c = size_ct(k)
print(f"Public Key Size: {total_bits_p} bits ({total_bytes_p} bytes)")
print(f"Ciphertext: {total_bits_c} bits ({total_bytes_c} bytes)")
