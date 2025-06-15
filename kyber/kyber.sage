#Class Kyber
class Kyber_KEM:
    def __init__(self, k_, n_, eta_1, eta_2, du_, dv_, q_):
        self.q = q_
        self.k = k_
        self.n = n_
        self.eta1 = eta_1
        self.eta2 = eta_2
        self.du = du_
        self.dv = dv_
    
    def keygen(self):
        A_ntt = create_matrix_ntt(self.n, self.k)
        s = create_vector(self.eta1, self.k)
        e = create_vector(self.eta1, self.k) 
        s_ntt = NTT_vector(s, self.n)
        e_ntt = NTT_vector(e, self.n)
        s_A_ntt = multiply_vector_matrix(s_ntt, A_ntt, self.k, self.n)
        t_ntt = sum_vectors(s_A_ntt, e_ntt, self.k) 
        return (A_ntt, t_ntt), s_ntt

    def encrypt(self, pk, m, rrho = 0):
        A, t = pk
        r = create_vector_(self.eta1, self.k, rrho)
        e1 = create_vector_(self.eta2, self.k, rrho) 
        e2 = create_vector_(self.eta2, 1, rrho)[0] 
        r_ntt = NTT_vector(r, self.n) 
        A_t = matrix_traspose(A, self.k)
        r_A = multiply_vector_matrix(r_ntt, A_t, self.k, self.n)
        r_A_intt = INTT_vector(r_A) 
        u = sum_vectors(r_A_intt, e1, self.k) 
        r_t_intt = INTT(multiply_vector_vector(r_ntt, t, self.k, self.n))
        v =  r_t_intt + e2 + (q // 2) * Rk(m)
        
        u_compressed = [mod_switching(x, 2^self.du, q) for x in u]
        v_compressed = mod_switching(v, 2^self.dv, q)

        return u_compressed, v_compressed

    def decrypt(self, sk, cipher):
        u, v = cipher  
        u_ = [mod_switching(x, q, 2^self.du) for x in u]
        v_ = mod_switching(v, q, 2^self.dv)
        u_ntt = NTT_vector(u_, self.n) 
        u_sk = multiply_vector_vector(u_ntt, sk, self.k, self.n)
        u_sk_intt = INTT(u_sk)
        m_dec = v_ - u_sk_intt
        m_dec_coeffs = mod_switching(m_dec, 2, self.q)

        return m_dec_coeffs
        
    def KEM_keygen(self):
        pk, sk = self.keygen()
        return pk, sk

    def KEM_encaps(self, pk):
        m = list(np.random.randint(0, 2, self.n))
        hash_ = H(m, pk)
        K, rrho = hash_[0:256], hash_[257:512]
        c = self.encrypt(pk, m, rrho)
        return hex(int(K, 2)), c

    def KEM_decaps(self, sk, c, K):
        m_prima = self.decrypt(sk, c)
        hash_ = H(m_prima, pk)
        K_prima, rrho_prima = hash_[0:256], hash_[257:512]
        c_prima = self.encrypt(pk, m_prima, rrho_prima)
        if c != c_prima:
            return None
        else:
            return hex(int(K_prima, 2))