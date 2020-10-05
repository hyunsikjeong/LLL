import itertools

class StateNotRecoveredException(Exception):
    pass

class CannotCalculateException(Exception):
    pass

class TLCGBreaker:

    def __init__(self, lcg_a, lcg_b, lcg_modulus_bit, num_truncated):
        self.lcg_a = lcg_a
        self.lcg_b = lcg_b
        # TODO: Support general modulus
        self.lcg_modulus_bit = lcg_modulus_bit
        self.lcg_modulus = 1 << lcg_modulus_bit
        self.state = None
        self.num_truncated = num_truncated

        self.output = []

    def set_outputs(self, outputs):
        self.output = outputs

    def get_output(self, idx):
        if self.state is None:
            raise StateNotRecoveredException()

        if idx < 0 and gcd(self.lcg_modulus, self.lcg_a) != 1:
            raise CannotCalculateException("Cannot calculate output with negative index" +\
                " when lcg_a and lcg_modulus is not coprime")

        if idx >= 0:
            state = self.state
            # TODO: O(lg idx)
            for _ in range(idx):
                state = (self.lcg_a * state + self.lcg_b) % self.lcg_modulus
        else:
            state = self.state
            inv_a = inverse_mod(self.lcg_a, self.lcg_modulus)
            inv_b = (-inv_a * self.lcg_b) % self.lcg_modulus
            # TODO: O(lg idx)
            for _ in range(-idx):
                state = (inv_a * state + inv_b) % self.lcg_modulus

        return int(state) >> self.num_truncated

    def recover_state(self):
        # 1. Initialize Lattice
        ln = len(self.output)
        L = Matrix(ZZ, ln - 1, ln - 1)

        for i in range(ln - 1):
            for j in range(ln - 1):
                if i == 0 and j == 0:
                    L[i, j] = self.lcg_modulus
                elif i == j:
                    L[i, j] = -1
                elif j == 0:
                    L[i, j] = pow(self.lcg_a, i, self.lcg_modulus)
                else:
                    L[i, j] = 0

        # 2. LLL!
        B = L.LLL()

        # 3. Find possible Xprime values
        # TODO: Make verbose setting to see the status
        for y_diff in itertools.product(range(2), repeat=ln - 1):
            Yprime = vector([ (self.output[i+1] - self.output[i] - y_diff[i]) & \
                ((1 << (self.lcg_modulus_bit - self.num_truncated)) - 1) for i in range(ln - 1)])
            v = vector([ round(RR(val) / self.lcg_modulus) * self.lcg_modulus - val \
                for val in 2^self.num_truncated * B * Yprime ])

            Zprime = B.solve_right(v)
            Xprime = 2^self.num_truncated * Yprime + Zprime

            # Check calculated Xprime is valid
            valid_Xprime = True
            for i in range(ln - 1):
                if int(Xprime[i]) >> self.num_truncated != Yprime[i]:
                    valid_Xprime = False
                    break

            g = gcd(self.lcg_a - 1, self.lcg_modulus)
            for i in range(ln - 1):
                if (Xprime[i] - self.lcg_b) % g != 0:
                    valid_Xprime = False
                    break

            if not valid_Xprime:
                continue

            print("YES")

            # Recover possible X from Xprime
            if self._recover_X(Xprime):
                return True

        # Failed to recover
        return False

    def _recover_X(self, Xprime):
        # Recover X from Xprime with binary search using
        # difference between real outputs and calculated outputs

        g = gcd(self.lcg_a - 1, self.lcg_modulus)
        x0 = (Xprime[0] - self.lcg_b) // g
        x0 = x0 * inverse_mod((self.lcg_a - 1) // g, self.lcg_modulus // g)
        x0 = x0 % (self.lcg_modulus // g)

        # TODO: Add option to get all possible states?
        # TODO: possible to add binary search using difference between outputs
        for i in range(g):
            x = x0 + (self.lcg_modulus // g) * i

            valid_x = True
            state = x
            for j in range(len(self.output)):
                if int(state) >> self.num_truncated != self.output[j]:
                    valid_x = False
                    break
                state = (self.lcg_a * state + self.lcg_b) % self.lcg_modulus

            if valid_x:
                self.state = x
                return True

        return False

# Example
if __name__ == "__main__":
    from Crypto.Util.number import bytes_to_long
    from os import urandom
    p = 2^32
    a = 0xdeadbeef
    b = 0x1337

    breaker = TLCGBreaker(a, b, 32, 20)

    state = bytes_to_long(urandom(4))
    output = []
    for i in range(100):
        prev = state
        state = (state * a + b) % p
        output.append(state >> 20)

    print(output)

    breaker.set_outputs(output[45:55])
    if breaker.recover_state():
        output = []
        for i in range(-45, 55):
            output.append(breaker.get_output(i))
        print(output)
    else:
        print("NO")
