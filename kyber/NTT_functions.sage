# SageMath implementation of polynomial multiplication using CRT in Z_q[X]/(X^d + α)
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
    
def NTT(poly, degree):
    list_ = []
    NTT_prev(poly, tree.root, degree, list_)
    return list_

def INTT(poly):
    result = INTT_prev(poly, tree.root, n)
    return result