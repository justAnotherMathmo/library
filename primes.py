from collections import defaultdict, deque
from random import randint

import numpy as np

from . import basic


class CreatePrimes(object):

    def __init__(self, bound, with_two=1):
        if with_two not in [0, 1]:
            raise ValueError('with2 must be boolean (1 or 0)')
        self.with_two = with_two

        # Gets all the primes up to bound as a np sieve
        sieve = np.ones(bound // 2, dtype=np.bool)
        for i in range(3, int(bound ** 0.5) + 1, 2):
            if sieve[i // 2]:
                sieve[i * i // 2::i] = False
        self.np_sieve =  2 * np.nonzero(sieve)[0][1::] + 1

    def as_array(self):
        if self.with_two:
            return np.hstack((2, self.np_sieve))
        else:
            return self.np_sieve

    def as_list(self):
        return self.as_array().tolist()

    def as_generator(self):
        return (p for p in self.as_array())

    def as_dict(self):
        result = defaultdict(bool)
        for val in self.as_array():
            result[val] = True
        return result

    def as_indexed_dict(self):
        result = {}
        index = 0
        for val in self.as_array():
            result[index] = int(val)
            index += 1
        return result

    def as_deque(self):
        return deque(self.as_array())


class PrimeCount(object):

    def __init__(self, m_bound=10**4, n_bound=10**6, **kwargs):
        # Some basic precomputation
        self.precomputed = {10 ** 11: 4118054813,
                            10 ** 12: 37607912018,
                            10 ** 13: 346065536839,
                            10 ** 14: 3204941750802,
                            10 ** 15: 29844570422669,
                            5 * 10 ** 15: 142377417196364,
                            10 ** 16: 279238341033925
                            }
        self.primes = kwargs.pop('primes', None)
        if self.primes is None:
            self.primes = CreatePrimes(n_bound, with_two=1)
        self.prime_array = self.primes.as_array()
        primes_gen = self.primes.as_generator()
        count = 0
        pNext = next(primes_gen)
        for n in range(n_bound + 1):
            if pNext == n:
                count += 1
                try:
                    pNext = next(primes_gen)
                except StopIteration:
                    pass
            self.precomputed[n] = count

        self.meissel_memoize = defaultdict(dict)
        for n in range(m_bound + 1):
            self.meissel_memoize[0][n] = n
        self.sieve = np.ones(m_bound + 1, dtype=np.bool)
        self.sieve[0] = False
        for i, p in enumerate(self.prime_array[:self.precomputed[m_bound]]):
            self.sieve[p::p] = False
            count = 0
            for n in range(1, m_bound + 1):
                if self.sieve[n]: count += 1
                self.meissel_memoize[i + 1][n] = count

        # Storing of variables that will be commonly used
        self.m_bound = m_bound
        self.n_bound = n_bound

    def __enter__(self):
        return self

    def __exit__(self):
        pass

    def _meissel_function_small(self, m, n):
        """The number of numbers <= m that are coprime to the first n prime numbers.
        Use for n < 1000 or so"""
        m = int(m)
        if n == 0:
            return m
        result = self.meissel_memoize[n].get(m, None)
        if result is None:
            result = self._meissel_function_small(m, n-1) - self._meissel_function_small(m / self.prime_array[n-1], n-1)
            self.meissel_memoize[n][m] = result
        return result

    def _meissel_function_large(self, m, n):
        """The number of numbers <= m that are coprime to the first n prime numbers.
        Run for larger values where repeating isn't going to happen often. Use for n > 1000 or so"""
        m = int(m)
        # if m <= 10000: return meis104[n][m]
        if n == 0:
            return m
        result = self.meissel_memoize[n].get(m, None)
        if result is None:
            result = 0
            primes_gen = (p for p in self.prime_array[:n:-1])
            stacks = defaultdict(lambda: defaultdict(int))
            stacks[n][m] = 1
            for N in range(n, 0, -1):
                prime_dividing = next(primes_gen)
                for M in stacks[N]:
                    # Marginal speed improvement?
                    # if M <= 10000 and n<1229:
                        # res += meis104[N][M]*stacks[N][M]
                        # continue
                    stacks[N-1][M] += stacks[N][M]
                    stacks[N-1][M//prime_dividing] -= stacks[N][M]
                del stacks[N]
            for M in stacks[0]:
                result += M*stacks[0][M]
        return result


    def primes_less_than(self, n):
        """Number of primes less than m using Lehmer Formula.
        Good for x<10^12 or so.
        Ensure that initPrimeCount() is run first."""
        n = int(n + 0.000000001)
        result = self.precomputed.get(n, None)
        if result is None:
            a = self.primes_less_than(pow(n, 1/4))
            c = self.primes_less_than(pow(n, 1/3))
            b = self.primes_less_than(pow(n, 1/2))
            result = self._meissel_function_small(n, a) + ((b + a - 2) * (b - a + 1)) // 2
            for i in range(a, b):
                result -= self.primes_less_than(n/self.prime_array[i])
                if i < c:
                    bi = self.primes_less_than(pow(n/self.prime_array[i], 1/2))
                    for j in range(i,bi):
                        result += j - self.primes_less_than(n/self.prime_array[i]/self.prime_array[j])
            self.precomputed[n] = result
        return result


class Factor(object):

    def __init__(self, primes):
        self.primes = primes

    def _f(self, x, N):
        """Psuedo-random number generator
        Helper function for _get_factor"""
        return (x**2 + 3) % N

    def small(self, num):
        """Prime factorises num (for small num)"""
        res = {}
        ma = num ** 0.5 + 1
        for i in self.primes.as_array():
            if i > ma:
                if num != 1: res[num] = 1
                break
            if num % i == 0:
                k = basic.pow_find(num, i)
                res[i] = k
                num = num // i ** k
                ma = num ** 0.5 + 1
        return res

    def _get_factor(self, N):
        """Algorithm for num < 2^70
        Algorithm described pages 7/8 in http://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf"""
        y, r, q, G = randint(N//4, 3*N//4), 1, 1, 1
        while G == 1:
            x = y
            for _ in range(r):
                y = self._f(y, N)
            k = 0
            while (G == 1) and (k < r):
                ys = y
                for _ in range(min(m, r-k)):
                    y = self._f(y, N)
                    q = (q * abs(x - y)) % N
                G = basic.gcd(q, N)
                k += m
            r *= 2
        if G == N:
            ys = y
            while True:
                ys = self._f(ys)
                G = basic.gcd(abs(x - ys), N)
                if G > 1:
                    break
        if G == N:
            raise UserWarning('Number not factored')
        return G



def is_prime_prob(num, accuracy=10):
    """Probabilistic prime checker (Miller-Rabin)"""
    s = basic.pow_find(num-1,2)
    d = (num-1)//2**s
    for kkk in range(accuracy):
        a = randint(num//4,3*num//4)
        if pow(a,d,num)==1: continue
        check = 1
        for r in range(s):
            if pow(a, 2**r*d, num) == num-1:
                check = 0
                break
        if check:
            return 0
    return 1


def legendre_symbol(n, p):
    """Gives the value (n/p). Note that p must be prime"""
    n = n%p
    if n==0: return 0
    if p==2: return 1
    if n==p-1:
        if p%4==1: return 1
        else: return -1
    if n==2:
        if p%8 in [1,7]: return 1
        else: return -1
    if round(n**0.5)**2==n: return 1
    factors = fact(n)
    res = 1
    for key in factors:
        if factors[key]%2==0: continue
        if key==2 and p%8 in [3,5]: res*=-1
        else:
            flippy = 1 - 2*(key%4==3)*(p%4==3)
            res *= legendre_symbol(p%key, key)*flippy
    return res


def tonelli_shanks(n, p):
    """Returns sqrt(n) mod p, with p a prime. Defined up to a minus sign"""
    L = legendre_symbol(n, p)
    if L == -1: raise ValueError("n has no solution or p isn't a prime")
    elif L == 0: return 0
    S = basic.pow_find(p-1, 2)
    Q = (p-1)//2**S
    if S == 1:
        return pow(n, (p+1)//4, p)
    z = 2
    while True:
        if legendre_symbol(z, p) == -1: break
        z = randint(3, p-1)
    c = pow(z, Q, p)
    R, t, M = pow(n, (Q+1)//2, p), pow(n, Q, p), S
    while True:
        if t == 1: return R
        tpow, i = pow(t, 2, p), 1
        while True:
            if tpow == 1: break
            i += 1
            tpow = pow(tpow, 2, p)
        b = pow(c, 2**(M-i-1), p)
        c = pow(b, 2, p)
        R, t, M = R*b%p, (t*c)%p, i