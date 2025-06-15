#KYBER
import numpy as np
import random
import time
import hashlib
import sys
import statistics

q, n = 3329, 256
R = IntegerModRing(q)  
PR.<X> = PolynomialRing(R)  
Rk = PR.quotient(X^n + 1)  
#Load Functions
load("class_tree.sage") 
load("NTT_functions.sage")
load("auxiliar_functions.sage")
load("kyber.sage")

#Load tree and inverse tree
tree, r_list = build_tree(1, n, q)
tree_inv = build_inverse_tree(tree, q, 1)

#Calculate size public key
def size_publickey(k_):
    bits_per_coef = ceil(log(q, 2))
    num_polynomials = k_^2 + k_
    total_coefficients = num_polynomials * n
    total_bits = total_coefficients * bits_per_coef
    total_bytes = ceil(total_bits / 8)
    return total_bits, total_bytes

#Calculate size ciphertext
def size_ct(k_):
    bits_per_coef = ceil(log(q, 2))
    num_polynomials = k_ + 1
    total_coefficients = num_polynomials * n
    total_bits = total_coefficients * bits_per_coef
    total_bytes = ceil(total_bits / 8)
    return total_bits, total_bytes

#Calculate size secret key
def size_secretkey(k_):
    bits_per_coef = ceil(log(q, 2))
    num_polynomials = 2 * k_ + 1   
    total_coefficients = num_polynomials * n
    total_bits = total_coefficients * bits_per_coef
    total_bytes = ceil(total_bits / 8)
    return total_bits, total_bytes

#Parameters Kyber
kyber_configs = [
    {"name": "Kyber-512", "k": 2, "eta1": 3, "eta2": 2, "du": 10, "dv": 4},
    {"name": "Kyber-768", "k": 3, "eta1": 2, "eta2": 2, "du": 10, "dv": 4},
    {"name": "Kyber-1024", "k": 4, "eta1": 2, "eta2": 2, "du": 11, "dv": 5},
]

#Testing time of every type of Kyber
for config in kyber_configs:
    print(f"############################ {config['name']} ############################")

    k = config["k"]
    #We create a kyber
    kyber = Kyber_KEM(k, n, config["eta1"], config["eta2"], config["du"], config["dv"], q)
    
    # Mesure time Keygen
    t0_keygen = time.process_time()
    pk, sk = kyber.KEM_keygen()
    t_keygen = time.process_time() - t0_keygen

    # Mesure time Encaps
    t0_encaps = time.process_time()
    K, c = kyber.KEM_encaps(pk)
    t_encaps = time.process_time() - t0_encaps

    # Mesure time Decaps 
    t0_decaps = time.process_time()
    k_shared = kyber.KEM_decaps(sk, c, K)
    t_decaps = time.process_time() - t0_decaps

    #Show results about time
    print("k_shared:", k_shared)
    print(f"KeyGen time   : {t_keygen:.6f} seconds")
    print(f"Encaps time   : {t_encaps:.6f} seconds")
    print(f"Decaps time   : {t_decaps:.6f} seconds")

    total_bits_p, total_bytes_p = size_publickey(k)
    total_bits_s, total_bytes_s = size_secretkey(k)
    total_bits_c, total_bytes_c = size_ct(k)

    #Print results about size
    print(f"Public Key Size : {total_bits_p} bits ({total_bytes_p} bytes)")
    print(f"Secret Key Size : {total_bits_s} bits ({total_bytes_s} bytes)")
    print(f"Ciphertext Size : {total_bits_c} bits ({total_bytes_c} bytes)\n")


num_runs = 10

results = []
#Official Benchmarking Time
for config in kyber_configs:
    keygen_times = []
    encaps_times = []
    decaps_times = []

    k = config["k"]
    #We create Kyber Class
    kyber = Kyber_KEM(k, n, config["eta1"], config["eta2"], config["du"], config["dv"], q)

    #Calculate statistics time of every version of Kyber
    for _ in range(num_runs):
        t0_keygen = time.process_time()
        pk, sk = kyber.KEM_keygen()
        t_keygen = time.process_time() - t0_keygen
        keygen_times.append(t_keygen)

        t0_encaps = time.process_time()
        K, c = kyber.KEM_encaps(pk)
        t_encaps = time.process_time() - t0_encaps
        encaps_times.append(t_encaps)

        t0_decaps = time.process_time()
        k_shared = kyber.KEM_decaps(sk, c, K)
        t_decaps = time.process_time() - t0_decaps
        decaps_times.append(t_decaps)

    keygen_avg = statistics.mean(keygen_times)
    keygen_med = statistics.median(keygen_times)
    encaps_avg = statistics.mean(encaps_times)
    encaps_med = statistics.median(encaps_times)
    decaps_avg = statistics.mean(decaps_times)
    decaps_med = statistics.median(decaps_times)

    results.append({
        "name": config["name"],
        "keygen_avg": keygen_avg,
        "keygen_med": keygen_med,
        "encaps_avg": encaps_avg,
        "encaps_med": encaps_med,
        "decaps_avg": decaps_avg,
        "decaps_med": decaps_med
    })

print("Version\tKeyGen Med (s)\tKeyGen Avg (s)\tEncaps Med (s)\tEncaps Avg (s)\tDecaps Med (s)\tDecaps Avg (s)")
for res in results:
    print(f"{res['name']}\t"
          f"{res['keygen_med']:.6f}\t{res['keygen_avg']:.6f}\t"
          f"{res['encaps_med']:.6f}\t{res['encaps_avg']:.6f}\t"
          f"{res['decaps_med']:.6f}\t{res['decaps_avg']:.6f}")

