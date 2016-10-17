import math
class Vec3D(object):
    def __init__(self, x, y,phi):
        self._x = x
        self._y = y
        self._phi = phi

    def __add__(self, other):
        return Vec3D(self._x + other._x, self._y + other._y,self._phi + other._phi)

    def __iadd__(self, other):
        if type(other) == int or type(other) == float:
            return Vec3D(self._x + other, self._y + other,self._phi + other)
        else:
            return Vec3D(self._x + other._x, self._y + other._y,self._phi + other._phi)

    def __sub__(self, other):
        return Vec3D(self._x - other._x, self._y - other._y,self._phi - other._phi)

    def __isub__(self, c):
        return Vec3D(self._x - c, self._y - c,self._phi - c)

    def __mul__(self, other):
        if type(other) == int or type(other) == float:
            return Vec3D(self._x*other, self._y*other, self._phi*other)
        else:
            return Vec3D(self._x*other._x, self._y*other._y, self._phi*other._phi)

    def __imul__(self, other):
        if type(other) == int or type(other) == float:
            return Vec3D(self._x*other, self._y*other, self._phi*other)
        else:
            return Vec3D(self._x*other._x, self._y*other._y, self._phi*other._phi)

    def __abs__(self):
        return math.sqrt(self._x**2 + self._y**2 + self._phi**2)

    def __eq__(self, other):
        return self._x == other._x and self._y == other._y and self._phi == other._phi

    def __str__(self):
        return '(%g, %g, %g)' % (self._x, self._y,self._phi)

    def __ne__(self, other):
        return not self.__eq__(other)  # reuse __eq__

#    @property
#    def x(self):
#        return self._x

#    @x.setter
#    def x(self,x):
#        self._x = x

#    @property
#    def y(self):
#        return self._y

#    @x.setter
#    def y(self,y):
#        self._y = y

#    @property
#    def phi(self):
#        return self._phi

#    @phi.setter
#    def phi(self,phi):
#        self._phi = phi
