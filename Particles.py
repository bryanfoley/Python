import graphics
from vector import Vec3D
from math import sqrt
colors = ['blue', 'orange', 'red', 'green', 'brown', 'purple', 'cyan', 'black']
G = Vec3D(0.0,-9.81,0.0)
class Particle(object):
    """Representation of a single particle in the simulation. A particle has a
    2D position, velocity, acceleration, frictional force and interacts with other nearby
    particles. Particles also maintain their graphical representation in the visualization."""
    density = 0.0005
    #cutoff = 0.01
    # prevent very large forces due to discretization/fp inaccuracy
    #min_r2 = (cutoff / 100) ** 2
    dt = 0.00005
    box_size = None
    scale_pos = None
    next_id = 0
    energy_correction = 1 # energy normalization

    def __init__(self, x, y,phi, vx, vy, vz, ax, ay, az,rad,mass,J,Z,mu,gamma,omega,Y,A):
        self.rtd0 = Vec3D(x,y,phi)
        self.rtd1 = Vec3D(vx,vy,vz)
        self.rtd2 = Vec3D(ax,ay,az)
        self.rtd3 = Vec3D(0.0,0.0,0.0)
        self.rtd4 = Vec3D(0.0,0.0,0.0)
        self.force = Vec3D(0.0,0.0,0.0)
        self.rad = rad
        self.mass = mass
        self.J = J
        self.Z = Z
        self.mu = mu
        self.gamma = gamma
        self.omega = omega
        self.Y = Y
        self.id = Particle.next_id
        Particle.next_id += 1
        self.graphic = None
        self.A = A

    def init_graphic(self, win, rad, owner=None):
        """Create a graphical representation of this particle for visualization.
        win is the graphics windown in which the particle should be drawn, rad
        is the radius of the particle, and owner is the thread/process number of
        the thread that owns this particle."""
        p = graphics.Point(self.rtd0._x * self.scale_pos + self.rad +5,
                           self.rtd0._y * self.scale_pos + self.rad +5)
        self.graphic = graphics.Circle(p, self.rad)
        color = colors[owner % len(colors)] if owner is not None else 'blue'
        self.graphic.setOutline(color)
        self.graphic.setFill(color)
        self.graphic.draw(win)

    def predict(self):
        self.oldx, self.oldy = self.rtd0._x, self.rtd0._y
        dt = self.dt
        a1 = dt
        a2 = a1 * dt/2
        a3 = a2 * dt/3
        a4 = a3 * dt/4

        self.rtd0 += ((self.rtd1*a1) + (self.rtd2*a2) + (self.rtd3*a3) + (self.rtd4*a4))
        self.rtd1 += self.rtd2*a2 + self.rtd3*a3 + self.rtd4*a4
        self.rtd2 += self.rtd3*a3 + self.rtd4*a4
        self.rtd3 += self.rtd4*a4

    def correct(self):
        self.corr = Vec3D(0.0,0.0,0.0)
        dt = self.dt 
        dt_recip = 1/dt

        coeff0 = (19.0)/(90.0)*(dt*dt/(2.0))
        coeff1 = (3.0)/(4.0)*(dt/(2.0));
        coeff3 = (1.0)/(2.0)*((3.0)*dt_recip)
        coeff4 = (1.0)/(12.0)*((12.0)*(dt_recip**2.0))

        self.accel=Vec3D(((1/self.mass)*self.force._x)+G._x,((1/self.mass)*self.force._y)+G._y,((1/self.J)*self.force._phi)+G._phi)
        self.corr = self.accel - self.rtd2

        self.rtd0 += self.corr*coeff0
        self.rtd1 += self.corr*coeff1
        self.rtd2 = self.accel
        self.rtd3 += self.corr*coeff3
        self.rtd4 += self.corr*coeff4

    def set_force_to_zero(self):
        self.force._x = 0.0
        self.force._y = 0.0
        self.force._phi = 0.0

    def apply_force(self,other):
        dx = other.rtd0._x - self.rtd0._x
        dy = other.rtd0._y - self.rtd0._y
        rr=sqrt((dx)**2 + (dy)**2);
        r1=self.rad;
        r2=other.rad;
        m1 = self.mass;
        m2 = other.mass;
        xi=(r1+r2)-rr;

        if(xi>0.0):
            self.Z+=1.0
            other.Z+=1.0
            #potential_sum += (1*(xi*xi))/2;
            Y=self.Y*other.Y/(self.Y+other.Y)
            A=0.5*(self.A+other.A)
            mu = self.mu if self.mu > other.mu else other.mu
            gamma = self.gamma if self.gamma > other.gamma else other.gamma
            gamma=1.0
            gamma_t = 1.0
            kn = 30.0
            reff = (self.rad*other.rad)/(self.rad+other.rad)
            meff = (self.mass*other.mass)/(self.mass+other.mass)
            dvx=self.rtd1._x-other.rtd1._x
            dvy=self.rtd1._y-other.rtd1._y
            rr_rez=1.0/rr
            ex=dx*rr_rez
            ey=dy*rr_rez
            xidot=-(ex*dvx+ey*dvy)
            vtrel=-dvx*ey + dvy*ex + self.omega*self.rad-other.omega*other.rad

            fn = sqrt(xi)*Y*sqrt(reff)*(xi+A*xidot)

            ft =- gamma_t*vtrel

            if fn < 0:
                fn=0.0
            if ft<-mu*fn:
                ft=-mu*fn
            if ft>mu*fn:
                ft=mu*fn

            self.force._x+=fn*ex-ft*ey
            self.force._y+=fn*ey+ft*ey
            self.force._phi+=self.rad*ft

            other.force._x+=-fn*ex+ft*ey
            other.force._y+=-fn*ey-ft*ey
            other.force._phi-=other.rad*ft

    def boundary(self):
        """Move a particle for one timestep. Slightly simplified Velocity Verlet
        integration conserves energy better than explicit Euler method."""
        #self.oldx, self.oldy = self.rtd0.x, self.rtd0.y
        #self.rtd1._x += self.rtd2._x * self.dt
        #self.rtd1._y += self.rtd2._y * self.dt
        #self.rtd0._x += self.rtd1._x * self.dt
        #self.rtd0._y += self.rtd1._y * self.dt

        # Bounce from walls
        size = self.box_size
        while self.rtd0._x < 0.0 or self.rtd0._x > size:
            self.rtd0._x = -self.rtd0._x if self.rtd0._x < 0.0 else 2 * size - self.rtd0._x
            self.rtd1._x = -self.rtd1._x
        while self.rtd0._y < 0.0 or self.rtd0._y > size:
            self.rtd0._y = -self.rtd0._y if self.rtd0._y < 0.0 else 2 * size - self.rtd0._y
            self.rtd1._y = -self.rtd1._y

    def move_graphic(self):
        """Move the assoicated graphic of this particle to its new location."""
        if self.graphic:
            dx, dy = self.rtd0._x - self.oldx, self.rtd0._y - self.oldy
            self.graphic.move(dx * self.scale_pos, dy * self.scale_pos)

    @property
    def energy(self):
        """Return the kinetic energy of this particle."""
        return 0.5 * self.mass * (self.rtd1._x ** 2 + self.rtd1._y ** 2)
    @property
    def penergy(self):
        """Return the potential energy of this particle"""
        pass

    def __repr__(self):
        fmt = "Particle({0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8})"
        return fmt.format(self.rtd0.x, self.rtd0.y, self.rtd0.phi, self.rtd1.x, self.rtd1.y, self.rtd1.phi, self.rtd2.x, self.rtd2.y, self.rtd2.phi)

    def __getstate__(self):
        """Remove graphic from state that is transferred to another process."""
        state = self.__dict__.copy()
        state['graphic'] = None
        return state
