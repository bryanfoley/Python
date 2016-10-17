from vector import Vec3D
from Particles import Particle
from random import random, seed, shuffle, gauss
from verlet import init_graphics
import unittest

class run_tests(unittest.TestCase):
    '''Vector class unit tests'''

    def test_001(self):
        '''A zero vector should be created'''
        vector1 = Vec3D(0.0,0.0,0.0)
        result = Vec3D(0.0,0.0,0.0)
        self.assertEqual(vector1,result)

    def test_002(self):
        '''A non-zero vector should be created'''
        vector1 = Vec3D(1.0,1.0,1.0)
        result = Vec3D(0.0,0.0,0.0)
        self.assertNotEqual(vector1,result)

    def test_003(self):
        '''Vector Addition'''
        vector1 = Vec3D(1.0,1.0,1.0)
        vector2 = Vec3D(2.0,2.0,2.0)
        result = Vec3D(3.0,3.0,3.0)
        vector3 = vector1+vector2
        self.assertEqual(vector3,result)

    def test_004(self):
        '''Vector Subtraction'''
        vector1 = Vec3D(1.0,1.0,1.0)
        vector2 = Vec3D(2.0,2.0,2.0)
        result = Vec3D(-1.0,-1.0,-1.0)
        vector3 = vector1-vector2
        self.assertEqual(vector3,result)

    def test_005(self):
        '''Vector Multiplication'''
        vector1 = Vec3D(1.0,1.0,1.0)
        vector2 = Vec3D(2.0,2.0,2.0)
        result = Vec3D(2.0,2.0,2.0)
        vector3 = vector1*vector2
        self.assertEqual(vector3,result)

    def test_006(self):
        '''Scalar Addition'''
        vector1 = Vec3D(1.0,1.0,1.0)
        scalar = 2.1
        result = Vec3D(3.1,3.1,3.1)
        vector1 += scalar
        self.assertEqual(vector1,result)

    def test_007(self):
        '''Scalar Subtraction'''
        vector1 = Vec3D(1.0,1.0,1.0)
        scalar = 2.1
        result = Vec3D(-1.1,-1.1,-1.1)
        vector1 -= scalar
        self.assertEqual(vector1,result)

    def test_008(self):
        '''Positive Scalar iMultiplication'''
        vector1 = Vec3D(1.0,1.0,1.0)
        scalar = 2.0
        result = Vec3D(2.0,2.0,2.0)
        vector1 *= scalar
        self.assertEqual(vector1,result)

    def test_009(self):
        '''Negative Scalar iMultiplication TC1'''
        vector1 = Vec3D(1.0,1.0,1.0)
        scalar = -3.0
        result = Vec3D(-3.0,-3.0,-3.0)
        vector1 *= scalar
        self.assertEqual(vector1,result)

    def test_010(self):
        '''Negative Scalar iMultiplication TC2'''
        vector1 = Vec3D(-1.0,1.0,1.0)
        scalar = -3.0
        result = Vec3D(3.0,-3.0,-3.0)
        vector1 *= scalar
        self.assertEqual(vector1,result)

    def test_011(self):
        '''Vector properties'''
        vector1 = Vec3D(-1.0,1.0,1.0)
        result = -1.0
        self.assertEqual(vector1._x,result)

    def test_012(self):
        '''Set Vector properties'''
        vector1 = Vec3D(-1.0,1.0,1.0)
        scalar1 = -2.0
        scalar2 = 5.1
        vector1._x = scalar1 + scalar2
        result = Vec3D(3.1,1.0,1.0)
        self.assertAlmostEqual(vector1, result, places=7, msg=None, delta=None)

    def test_013(self):
        '''Augment Vector properties'''
        vector1 = Vec3D(-1.0,1.0,1.0)
        scalar1 = -2.0
        scalar2 = 5.1
        vector1._x += scalar1 + scalar2
        result = Vec3D(2.1,1.0,1.0)
        self.assertAlmostEqual(vector1, result, places=7, msg=None, delta=None)

    def test_014(self):
        '''Augment Vector'''
        vector1 = Vec3D(-1.0,1.0,1.0)
        scalar1 = -2.0
        vector2 = Vec3D(-2.0,0.0,5.1)
        vector1 += vector2*scalar1
        result = Vec3D(3.0,1.0,-9.2)
        self.assertAlmostEqual(vector1, result, places=7, msg=None, delta=None)

    '''Particle class unit tests'''

if __name__ == '__main__':
    unittest.main()
