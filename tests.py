# Python Libraries
import types
import unittest
from collections import deque
from math import cos, pi, sin

# External Libraries
import numpy as np

# Custom Libraries
from . import primes
from . import basic


class TestBasicFunctions(unittest.TestCase):

    def test_gcd(self):
        self.assertEqual(
            basic.gcd(3, 4),
            1
        )
        self.assertEqual(
            basic.gcd(6, 10),
            2
        )
        self.assertEqual(
            basic.gcd(80, 0),
            80
        )
        self.assertEqual(
            basic.gcd(6, 22),
            2
        )
        self.assertEqual(
            basic.gcd(2**5*3**24*5**8, 2**20*3**13*5**7),
            2**5*3**13*5**7
        )

    def test_gcd_coeff(self):
        x, y = basic.gcd_coeff(4, 5)
        self.assertEqual(
            x*4 + y*5,
            1
        )
        x, y = basic.gcd_coeff(32, 48)
        self.assertEqual(
            x*32 + y*48,
            16
        )
        x, y = basic.gcd_coeff(146714362, 24352627)
        self.assertEqual(
            x*146714362 + y*24352627,
            1
        )

    def test_mod_inv(self):
        self.assertEqual(
            basic.mod_inv(1, 5),
            1
        )
        self.assertEqual(
            basic.mod_inv(3, 4),
            3
        )
        self.assertEqual(
            basic.mod_inv(4, 13),
            10
        )
        self.assertEqual(
            basic.mod_inv(57, 149),
            34
        )

    def test_pow_find(self):
        self.assertEqual(
            basic.pow_find(2, 2),
            1
        )
        self.assertEqual(
            basic.pow_find(48, 2),
            4
        )
        self.assertEqual(
            basic.pow_find(2**4*3**7*5**12*7**9, 6),
            4
        )
        self.assertEqual(
            basic.pow_find(1, 6),
            0
        )
        with self.assertRaises(ValueError):
            basic.pow_find(0, 2)

    def test_int_sqrt(self):
        self.assertEqual(
            basic.int_sqrt(16),
            4
        )
        self.assertEqual(
            basic.int_sqrt(24),
            4
        )
        self.assertEqual(
            basic.int_sqrt(24, floor=0),
            5
        )
        self.assertEqual(
            basic.int_sqrt(6332850846508003265186595),
            2516515616185
        )
        self.assertEqual(
            basic.int_sqrt(6332850846508003265186596),
            2516515616186
        )

    def test_newton_rhapson(self):
        self.assertAlmostEqual(
            basic.newton_rhapson(cos, lambda x: -sin(x), 1),
            pi/2,
            places=9
        )

    def test_binary_search(self):
        self.assertAlmostEqual(
            basic.binary_search(cos,0, pi),
            pi/2,
            places=9
        )

class TestCreatePrimes(unittest.TestCase):

    def setUp(self):
        self.primes = primes.CreatePrimes(30)
        self.primes_no_two = primes.CreatePrimes(30, with_two=0)
        self.correct_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        self.correct_list_no_two = [3, 5, 7, 11, 13, 17, 19, 23, 29]

    def test_array(self):
        self.assertIsInstance(
            self.primes.as_array(),
            np.ndarray
        )
        np.testing.assert_array_equal(
            self.primes.as_array(),
            np.array(self.correct_list)
        )
        np.testing.assert_array_equal(
            self.primes_no_two.as_array(),
            np.array(self.correct_list_no_two)
        )

    def test_list(self):
        self.assertIsInstance(
            self.primes.as_list(),
            list
        )

    def test_generator(self):
        self.assertIsInstance(
            self.primes.as_generator(),
            types.GeneratorType
        )

    def test_dict(self):
        dictionary = self.primes.as_dict()
        self.assertIsInstance(
            dictionary,
            dict
        )
        for key, value in [(0, 0), (2, 1), (3, 1), (6, 0), (11, 1), (10**50, 0)]:
            self.assertEqual(
                dictionary[key],
                value
            )
        self.assertEqual(
            self.primes_no_two.as_dict()[2],
            0
        )


    def test_deque(self):
        self.assertIsInstance(
            self.primes.as_deque(),
            deque
        )

class TestFactor(unittest.TestCase):

    def setUp(self):
        self.primes = primes.CreatePrimes()
        self.factor = primes.Factor()

    def test__f(self):
        self.assertEqual(
            self.factor._f(4, 7),
            5
        )

    def test_small(self):
        pass

    def test__get_factor(self):
        pass


class TestPrimeCount(unittest.TestCase):

    def setUp(self):
        self.Counter = primes.PrimeCount()

    def test_meissel_function(self):
        for meissel_function in [self.Counter._meissel_function_small,
                                 self.Counter._meissel_function_large]:
            self.assertEqual(
                meissel_function(10, 2),
                3
            )
            self.assertEqual(
                meissel_function(50, 2),
                17
            )
            self.assertEqual(
                meissel_function(50, 3),
                14
            )
            self.assertEqual(
                meissel_function(10**5, 40),
                10287
            )

    def test_primes_less_than(self):
        self.assertEqual(
            self.Counter.primes_less_than(100),
            25
        )
        self.assertEqual(
            self.Counter.primes_less_than(726445883),
            37550272
        )
        self.assertEqual(
            self.Counter.primes_less_than(123456),
            11601
        )

class TestMiscPrimeFunctions(unittest.TestCase):

    def test_is_prime_prob(self):
        """This has a 1 in 2^40 ~ 10^12 chance of failing"""
        self.assertTrue(
            primes.is_prime_prob(2011, accuracy=20)
        )
        self.assertTrue(
            primes.is_prime_prob(899809363, accuracy=20)
        )
        self.assertFalse(
            primes.is_prime_prob(7081, accuracy=20)
        )
        with self.assertRaises(ValueError):
            primes.is_prime_prob(1, accuracy=20)

    def test_legendre_symbol(self):
        self.assertEqual(
            primes.legendre_symbol(16, 53),
            1
        )
        self.assertEqual(
            primes.legendre_symbol(21, 29),
            -1
        )
        self.assertEqual(
            primes.legendre_symbol(7, 17),
            -1
        )
        self.assertEqual(
            primes.legendre_symbol(15, 127),
            1
        )
        self.assertEqual(
            primes.legendre_symbol(29+5000*83, 83),
            1
        )
        self.assertEqual(
            primes.legendre_symbol(146, 73),
            0
        )


    def test_tonelli_shanks(self):
        self.assertIn(
            primes.tonelli_shanks(25, 73),
            [5, 68]
        )
        self.assertIn(
            primes.tonelli_shanks(1273, 2011),
            [861, 1150]
        )
        self.assertIn(
            primes.tonelli_shanks(656518566, 899809363),
            [12345678, 887463685]
        )
        with self.assertRaises(ValueError):
            primes.tonelli_shanks(7, 17)


if __name__ == '__main__':
    unittest.main()