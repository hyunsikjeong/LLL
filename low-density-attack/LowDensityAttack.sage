class HighDensityException(Exception):
    pass


class LOAttack:

    def __init__(self, array, target_sum, try_on_high_density=False):
        self.array = array
        self.n = len(self.array)
        self.target_sum = target_sum
        self.density = self._calc_density()
        self.try_on_high_density = try_on_high_density

    def _calc_density(self):
        return self.n / log(max(self.array), 2)

    def _check_ans(self, ans):
        calc_sum = sum(map(lambda x: x[0] * x[1], zip(self.array, ans)))
        return self.target_sum == calc_sum

    def solve(self):
        if self.density >= 0.6463 and not self.try_on_high_density:
            raise HighDensityException()

        # 1. Initialize Lattice
        L = Matrix(ZZ, self.n + 1, self.n + 1)
        N = ceil(self.n ^ 0.5 / 2)
        for i in range(self.n + 1):
            for j in range(self.n + 1):
                if j == self.n and i < self.n:
                    L[i, j] = N * self.array[i]
                elif j == self.n:
                    L[i, j] = N * self.target_sum
                elif i == j:
                    L[i, j] = 1
                else:
                    L[i, j] = 0

        # 2. LLL!
        B = L.LLL()

        # 3. Find answer
        for i in range(self.n + 1):
            if B[i, self.n] != 0:
                continue

            if all(-1 <= v <= 0 for v in B[i]):
                ans = [-B[i, j] for j in range(self.n)]
                if self._check_ans(ans):
                    return ans

        # Failed to find answer
        return None


class CJLOSSAttack:

    def __init__(self, array, target_sum, try_on_high_density=False):
        self.array = array
        self.n = len(self.array)
        self.target_sum = target_sum
        self.density = self._calc_density()
        self.try_on_high_density = try_on_high_density

    def _calc_density(self):
        return self.n / log(max(self.array), 2)

    def _check_ans(self, ans):
        calc_sum = sum(map(lambda x: x[0] * x[1], zip(self.array, ans)))
        return self.target_sum == calc_sum

    def solve(self):
        if self.density >= 0.9408 and not self.try_on_high_density:
            raise HighDensityException()

        # 1. Initialize Lattice
        L = Matrix(ZZ, self.n + 1, self.n + 1)
        N = ceil(self.n ^ 0.5 / 2)
        for i in range(self.n + 1):
            for j in range(self.n + 1):
                if j == self.n and i < self.n:
                    L[i, j] = 2 * N * self.array[i]
                elif j == self.n:
                    L[i, j] = 2 * N * self.target_sum
                elif i == j:
                    L[i, j] = 2
                elif i == self.n:
                    L[i, j] = 1
                else:
                    L[i, j] = 0

        # 2. LLL!
        B = L.LLL()

        # 3. Find answer
        for i in range(self.n + 1):
            if B[i, self.n] != 0:
                continue

            if all(v == -1 or v == 1 for v in B[i][:self.n]):
                ans = [ (-B[i, j] + 1) // 2 for j in range(self.n)]
                if self._check_ans(ans):
                    return ans

        # Failed to find answer
        return None


# Example
if __name__ == "__main__":
    from sage.misc.prandom import randint

    def gen_example(n, density):
        max_val = floor(2 ^ (n / density))
        arr = [randint(1, max_val) for _ in range(n)]
        ans = [randint(0, 1) for _ in range(n)]
        target_sum = sum(map(lambda x: x[0] * x[1], zip(arr, ans)))

        return arr, target_sum, ans

    # 1. Density == 0.5 (Solvable by both)
    print("Density: 0.5")
    arr, target_sum, ans = gen_example(32, 0.5)
    print("ans:   ", ans)
    attack = LOAttack(arr, target_sum, True)
    print("LO:    ", attack.solve())
    attack = CJLOSSAttack(arr, target_sum, True)
    print("CJLOSS:", attack.solve())

    # 2. Density == 0.7 (Only solvable by CJLOSS)
    print("Density: 0.7")
    arr, target_sum, ans = gen_example(32, 0.7)
    print("ans:   ", ans)
    attack = LOAttack(arr, target_sum, True)
    print("LO:    ", attack.solve())
    attack = CJLOSSAttack(arr, target_sum, True)
    print("CJLOSS:", attack.solve())

    # 3. Density == 0.9 (Only solvable by CJLOSS)
    print("Density: 0.9")
    arr, target_sum, ans = gen_example(32, 0.9)
    print("ans:   ", ans)
    attack = LOAttack(arr, target_sum, True)
    print("LO:    ", attack.solve())
    attack = CJLOSSAttack(arr, target_sum, True)
    print("CJLOSS:", attack.solve())

    # 4. Density == 0.95 (Not solvable by both)
    print("Density: 0.95")
    arr, target_sum, ans = gen_example(32, 0.95)
    print("ans:   ", ans)
    attack = LOAttack(arr, target_sum, True)
    print("LO:    ", attack.solve())
    attack = CJLOSSAttack(arr, target_sum, True)
    print("CJLOSS:", attack.solve())
