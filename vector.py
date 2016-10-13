class Vec3D(object):
    def __init__(self, x, y,phi):
        self._x = x
        self._y = y
        self._phi = phi

    def __add__(self, other):
        return Vec3D(self._x + other.x, self._y + other.y,self._phi + other.phi)

    def __sub__(self, other):
        return Vec3D(self._x - other.x, self._y - other.y,self._phi - other.phi)

    def __mul__(self, other):
        return self._x*other.x + self._y*other.y + self._phi*other.phi

    def __mul__(self, c):
        return self._x*c + self._y*c + self._phi*c

    def __abs__(self):
        return math.sqrt(self._x**2 + self._y**2 + self._phi**2)

    def __eq__(self, other):
        return self._x == other.x and self._y == other.y and self._phi == other.phi

    def __str__(self):
        return '(%g, %g %g)' % (self._x, self._y,self._phi)

    def __ne__(self, other):
        return not self.__eq__(other)  # reuse __eq__

    def setx(self, x):
        self._x = x

    def sety(self, y):
        self._y = y

    def setphi(self, phi):
        self._phi = phi

    def setvec(self,other):
        self._x = other.x
        self._y = other.y
        self._phi = other.phi

    def getx(self):
        return self._x

    def gety(self):
        return self._y

    def getphi(self):
        return self._phi
