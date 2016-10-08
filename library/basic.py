def gcd(a, b):
    """Greatest common divisor of a and b"""
    while a != 0:
        a, b = b%a, a
    return b


def gcd_coeff(a, b):
    """Gives the coefficients such that x*a + y*b = gcd(a, b)"""
    u1, u2 = 1, 0
    v1, v2 = 0, 1
    while a != 0:
        q = b // a
        a, b = b % a, a
        u1, u2 = u2, u1 - q * u2
        v1, v2 = v2, v1 - q * v2
    return v1, u1


def mod_inv(num, mod):
    """Finds num^-1 modulo mod"""
    r1, r2, a, b= num%mod, mod, 1, 0
    while r1!=1:
        r1, r2, a, b = r2%r1, r1, b-a*(r2//r1), a
    return a%mod


def pow_find(num, div):
    """Finds s such that div^s||num"""
    if num == 0: raise ValueError('num must be a positive integer')
    if num % div != 0: return 0
    count = 0
    power_check = div
    while True:
        if num % power_check != 0: break
        power_check *= div
        count += 1
    return count


def int_sqrt(n, floor = 1):
    """Largest integer smaller than sqrt(n), ceiling function if floor = 0"""
    x = n
    y = (x + 1) >> 1
    while y < x:
        x = y
        y = (x + n // x) >> 1
    if floor: return x
    else:
        if x**2 < n:
            return x+1
        return x


def newton_rhapson(func, func_d, x0, tol=10**(-9), itera=None):
    """Uses the Newton Rhapson method to find an approximate root to func near x0.
    func_d is the derivative of func.
    If itera is not None, will limit to itera iterations, as well as tolerance"""
    res = x0
    loops = 0
    while True:
        res = res - func(res)/func_d(res)
        loops += 1
        if abs(func(res)) < tol or (itera is not None and loops >= itera):
            break
    return res


def binary_search(func, x_min, x_max, tol=10**(-13), itera=None):
    """Finds an approximate root to func within the range x_min to x_max via a binary search
    If itera is not None, will limit to itera iterations, as well as tolerance"""
    l_sign, r_sign = 2*(func(x_min) > 0) - 1, 2*(func(x_max) > 0) -1
    if l_sign == r_sign:
        raise ValueError("func(x_min) and func(x_max) must have opposite signs")
    while abs(func(x_max)) > tol:
        x_new = (x_min + x_max)/2
        n_sign = 2*(func(x_new) > 0) - 1
        if n_sign == l_sign:
            x_min = x_new
        else:
            x_max = x_new
    return x_new
