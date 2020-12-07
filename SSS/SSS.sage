# Code base from https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/boneh_durfee.sage
# Based on http://souravsengupta.com/publications/2010_indocrypt_2.pdf

import time

############################################
# Config
##########################################

"""
Setting debug to true will display more informations
about the lattice, the bounds, the vectors...
"""
debug = True

"""
Setting strict to true will stop the algorithm (and
return (-1, -1)) if we don't have a correct 
upperbound on the determinant. Note that this 
doesn't necesseraly mean that no solutions 
will be found since the theoretical upperbound is
usualy far away from actual results. That is why
you should probably use `strict = False`
"""
strict = False

"""
This is experimental, but has provided remarkable results
so far. It tries to reduce the lattice as much as it can
while keeping its efficiency. I see no reason not to use
this option, but if things don't work, you should try
disabling it
"""
helpful_only = True
dimension_min = 7 # stop removing if lattice reaches that dimension

############################################
# Functions
##########################################

# display stats on helpful vectors
def helpful_vectors(BB, modulus):
    nothelpful = 0
    for ii in range(BB.dimensions()[0]):
        if BB[ii,ii] >= modulus:
            nothelpful += 1

    print(nothelpful, "/", BB.dimensions()[0], " vectors are not helpful")

# display matrix picture with 0 and X
def matrix_overview(BB, bound):
    for ii in range(BB.dimensions()[0]):
        a = ('%02d ' % ii)
        for jj in range(BB.dimensions()[1]):
            a += '0' if BB[ii,jj] == 0 else 'X'
            if BB.dimensions()[0] < 60:
                a += ' '
        if BB[ii, ii] >= bound:
            a += '~'
        print(a)

# tries to remove unhelpful vectors
# we start at current = n-1 (last vector)
def remove_unhelpful(BB, monomials, bound, current):
    # end of our recursive function
    if current == -1 or BB.dimensions()[0] <= dimension_min:
        return BB

    # we start by checking from the end
    for ii in range(current, -1, -1):
        # if it is unhelpful:
        if BB[ii, ii] >= bound:
            affected_vectors = 0
            affected_vector_index = 0
            # let's check if it affects other vectors
            for jj in range(ii + 1, BB.dimensions()[0]):
                # if another vector is affected:
                # we increase the count
                if BB[jj, ii] != 0:
                    affected_vectors += 1
                    affected_vector_index = jj

            # level:0
            # if no other vectors end up affected
            # we remove it
            if affected_vectors == 0:
                print("* removing unhelpful vector", ii)
                BB = BB.delete_columns([ii])
                BB = BB.delete_rows([ii])
                monomials.pop(ii)
                BB = remove_unhelpful(BB, monomials, bound, ii-1)
                return BB

            # level:1
            # if just one was affected we check
            # if it is affecting someone else
            elif affected_vectors == 1:
                affected_deeper = True
                for kk in range(affected_vector_index + 1, BB.dimensions()[0]):
                    # if it is affecting even one vector
                    # we give up on this one
                    if BB[kk, affected_vector_index] != 0:
                        affected_deeper = False
                # remove both it if no other vector was affected and
                # this helpful vector is not helpful enough
                # compared to our unhelpful one
                if affected_deeper and abs(bound - BB[affected_vector_index, affected_vector_index]) < abs(bound - BB[ii, ii]):
                    print("* removing unhelpful vectors", ii, "and", affected_vector_index)
                    BB = BB.delete_columns([affected_vector_index, ii])
                    BB = BB.delete_rows([affected_vector_index, ii])
                    monomials.pop(affected_vector_index)
                    monomials.pop(ii)
                    BB = remove_unhelpful(BB, monomials, bound, ii-1)
                    return BB
    # nothing happened
    return BB


""" 
Returns:
* 0,0   if it fails
* -1,-1 if `strict=true`, and determinant doesn't bound
* x0,y0 the solutions of `pol`
"""
def attack(N, l0, e, m, t, X, Y):
    modulus = e

    PR.<x, y> = PolynomialRing(ZZ)
    fMSB = 1 + (x + l0) * (N + y)

    # x-shifts
    gg = []
    for j in range(m + 1):
        for i in range(m - j + t + 1):
            xshift = x^i * fMSB^j * e^(m - j)
            gg.append(xshift)

    # y-shifts    
    for j in range(m + 1):
        for i in range(1, m - j + 1):
            yshift = y^i * fMSB^j * e^(m - j)
            gg.append(yshift)

    def order_gg(idx, gg, monomials):
        if idx == len(gg):
            return gg, monomials

        for i in range(idx, len(gg)):
            polynomial = gg[i]
            non = []
            for monomial in polynomial.monomials():
                if monomial not in monomials:
                    non.append(monomial)
            
            if len(non) == 1:
                new_gg = gg[:]
                new_gg[i], new_gg[idx] = new_gg[idx], new_gg[i]

                return order_gg(idx + 1, new_gg, monomials + non)    

    gg, monomials = order_gg(0, gg, [])

    # construct lattice B
    nn = len(monomials)
    BB = Matrix(ZZ, nn)
    for ii in range(nn):
        BB[ii, 0] = gg[ii](0, 0)
        for jj in range(1, nn):
            if monomials[jj] in gg[ii].monomials():
                BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj]) * monomials[jj](X, Y)

    # Prototype to reduce the lattice
    if helpful_only:
        # automatically remove
        BB = remove_unhelpful(BB, monomials, modulus^m, nn-1)
        # reset dimension
        nn = BB.dimensions()[0]
        if nn == 0:
            print("failure")
            return 0,0

    # check if vectors are helpful
    if debug:
        helpful_vectors(BB, modulus^m)
    
    # check if determinant is correctly bounded
    det = BB.det()
    bound = modulus^(m*nn)
    if det >= bound:
        print("We do not have det < bound. Solutions might not be found.")
        print("Try with highers m and t.")
        if debug:
            diff = (log(det) - log(bound)) / log(2)
            print("size det(L) - size e^(m*n) = ", floor(diff))
        if strict:
            return -1, -1
    else:
        print("det(L) < e^(m*n) (good! If a solution exists < N^delta, it will be found)")

    # display the lattice basis
    if debug:
        matrix_overview(BB, modulus^m)

    # LLL
    if debug:
        print("optimizing basis of the lattice via LLL, this can take a long time")

    BB = BB.LLL()

    if debug:
        print("LLL is done!")

    # transform vector i & j -> polynomials 1 & 2
    if debug:
        print("looking for independent vectors in the lattice")
    found_polynomials = False
    
    for pol1_idx in range(nn - 1):
        for pol2_idx in range(pol1_idx + 1, nn):
            # for i and j, create the two polynomials
            PR.<a, b> = PolynomialRing(ZZ)
            pol1 = pol2 = 0
            for jj in range(nn):
                pol1 += monomials[jj](a,b) * BB[pol1_idx, jj] / monomials[jj](X, Y)
                pol2 += monomials[jj](a,b) * BB[pol2_idx, jj] / monomials[jj](X, Y)

            # resultant
            PR.<q> = PolynomialRing(ZZ)
            rr = pol1.resultant(pol2)

            # are these good polynomials?
            if rr.is_zero() or rr.monomials() == [1]:
                continue
            else:
                print("found them, using vectors", pol1_idx, "and", pol2_idx)
                found_polynomials = True
                break
        if found_polynomials:
            break

    if not found_polynomials:
        print("no independant vectors could be found. This should very rarely happen...")
        return 0, 0
    
    rr = rr(q, q)

    # solutions
    soly = rr.roots()

    if len(soly) == 0:
        print("Your prediction (delta) is too small")
        return 0, 0

    soly = soly[0][0]
    ss = pol1(q, soly)
    solx = ss.roots()[0][0]

    return solx, soly

def example(beta, delta, m=8, t=0):
    from sage.misc.prandom import randint

    _lambda = max(delta, beta - 1/2)

    tmp = 1/12 * m^3 - 13/12 * m + 1/4 * m^2 * t + 1/4 * m * t
    tmp /= 1/2 * m^3 + m^2 + 1/2 * m + 1/2 * t^2 + 1/2 * t + m^2 * t + 1/2 * m * t^2 + 3/2 * m * t

    if _lambda >= tmp:
        print(f"assertion fail on beta: {beta}, delta: {delta}, m: {m}, t: {t}")
        print("Try with another m & t")
        return False

    nbits = 1024

    p = random_prime(2^(nbits // 2))
    q = random_prime(2^(nbits // 2))

    N = p * q
    phi = (p - 1) * (q - 1)

    while True:
        d = randint(N ** (beta - 0.001), N ** beta)
        if gcd(d, phi) == 1:
            break

    e = inverse_mod(d, phi)

    d1_bits = floor(delta * nbits)
    M = 2 ^ d1_bits
    d1, d0 = d % M, d // M

    print("beta =", beta)
    print("delta =", delta)
    print("N =", N)
    print("d =", d)
    print("e =", e)
    print("M =", M)
    print("d0 =", d0)
    print("d1 =", d1)

    X = floor(N^delta)
    Y = floor(N^(1/2))

    l0 = floor((e * d0 * M - 1) / N)

    solx, soly = attack(N, l0, e, m, t, X, Y)

    # found a solution?
    if solx > 0:
        print("=== solution found ===")
        if False:
            print("x:", solx)
            print("y:", soly)

        tmp = 1 + (solx + l0) * (N + soly)
        d = int(tmp / e)
        print("private key found:", d)
        return True
    else:
        print("=== no solution was found ===")
        return False

if __name__ == "__main__":

    arr = [
        (0.3, 0.15)
    ]

    for v1, v2 in arr:
        found = False
        for m in range(5, 20):
            for t in range(m):
                if example(v1, v2, m, t):
                    print("NICE!!!!!!!!!!!!!!!!")
                    found = True
                    break
            
            if found:
                break
        
        print("=====================================================")