from random import random, seed, shuffle, gauss
import graphics
from math import ceil, sqrt
from functools import reduce
from time import time
from ucb import main
from Particles import Particle

default_num_particles = 20
default_steps = 10000

##################
# Initialization #
##################

def make_particles(n):
    """Construct a list of n particles in two dimensions, initially distributed
    evenly but with random velocities. The resulting list is not spatially
    sorted."""
    seed(1000)
    sx = ceil(sqrt(n))
    sy = (n + sx - 1) // sx
    start_id = Particle.next_id
    Particle.box_size = sqrt(Particle.density * n)
    particles = [Particle(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0,0.0,0.0,gauss(10.5,0.5),1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0) for _ in range(n)]
    size = Particle.box_size

    # Make sure particles are not spatially sorted
    shuffle(particles)

    for p in particles:
        # Distribute particles evenly to ensure proper spacing
        i = p.id - start_id
        p.rtd0._x = size * (1 + i % sx) / (1 + sx)
        p.rtd0._y = size * (1 + i / sx) / (1 + sy)

        # Assign random velocities within a bound
        p.rtd1._x = random() * 2 - 1
        p.rtd1._y = random() * 2 - 1

    return particles

def init_graphics(distribution, total, update_interval=1, size=600):
    """Initialize the visualization, if update_interval is nonzero. distribution
    is the set of particles, divided into lists for each thread or process.
    total is the total number of particles. size is the base size of the
    simulation; the window size will be slightly larger."""
    if not update_interval:
        return None, None

    psize = ceil(sqrt(10000 / total)) # particle size
    # Adjust window size so that particle edges do not go off the screen
    wsize = size + psize * 2 + 5
    win = graphics.GraphWin('Particle Simulation', wsize, wsize,
                            autoflush=False)
    win.setBackground('white')

    # Initialize particle graphics
    Particle.scale_pos = size / Particle.box_size
    energy = 0
    for t in range(len(distribution)):
        particles = distribution[t]
        for p in particles:
            p.init_graphic(win, psize, t)
            energy += p.energy

    # Initialize step number
    text = graphics.Text(graphics.Point(wsize // 2, 20),
                         'step = 0, energy = ' + str(energy))
    text.setSize(18)
    text.draw(win)

    return win, text

def update_step(win, text, step, energy, update_interval):
    """Update the visualization if appropriate given the step number and update
    interval."""
    if update_interval and step % update_interval == 0:
        format_str = 'step = {0}, energy = {1}'
        text.setText(format_str.format(step, round(1000 * energy)))
        win.update()

#####################
# Serial Simulation #
#####################

def serial_simulation(n, steps, num_threads=1, normalize_energy=False,
                      update_interval=1):
    """Simulate n particles sequentially for steps steps. num_threads should
    always be 1. update_interval is the visualization update interval."""
    assert num_threads == 1, 'serial_simulation cannot use multiple threads'

    # Create particles
    particles = make_particles(n)
    #initial_energy = reduce(lambda x, p: x + p.energy, particles, 0)

    # Initialize visualization
    win, text = init_graphics((particles,), n, update_interval)

    # Perform simulation
    start = time()
    for step in range(steps):
        # Compute forces
        for p1 in particles:
           # p1.rtd2._x = p1.rtd2._y = 0 # reset accleration to 0
            p1.set_force_to_zero()
            p1.predict()
            p1.boundary()
            for p2 in particles:
                if p2.id is not p1.id:
                    p1.apply_force(p2)
            p1.correct()

        # Move particles
        for p in particles:
            # Energy normalization
            p.rtd1._x *= Particle.energy_correction
            p.rtd1._y *= Particle.energy_correction

        # Update visualization
        energy = 0
        for p in particles:
            p.move_graphic()
            energy += p.energy
        update_step(win, text, step, energy, update_interval)

        # Energy normalization
        if normalize_energy:
            Particle.energy_correction = sqrt(initial_energy / energy)
    end = time()

    print('serial simulation took {0} seconds'.format(end - start))

@main
def run(*args):
    simulation, num_threads = serial_simulation, 1
    num_particles, steps = default_num_particles, default_steps
    normalize_energy = False
    update_interval = 0
    i = 0
    while i < len(args):
        if args[i] == '-t':
            simulation = thread_simulation
            num_threads = int(args[i+1])
        elif args[i] == '-p':
            simulation = process_simulation
            num_threads = int(args[i+1])
        elif args[i] == '-n':
            num_particles = int(args[i+1])
        elif args[i] == '-s':
            steps = int(args[i+1])
        elif args[i] == '-g' or args[i] == '-v':
            update_interval = 1
            i -= 1
        elif args[i] == '-u':
            update_interval = int(args[i+1])
        elif args[i] == '-dt':
            Particle.dt = float(args[i+1])
        elif args[i] == '-e':
            normalize_energy = True
            i -= 1
        else:
            if args[i] != '-h' and args[i] != '-help':
                print('unknown argument:', args[i], file=sys.stderr)
            print('Options:\n' +
                  '  -t <num>     run with <num> threads\n' +
                  '  -p <num>     run with <num> processes\n' +
                  '  -n <num>     simulate <num> particles\n' +
                  '  -s <num>     run for <num> timesteps\n' +
                  '  -v, -g       enable visualization\n' +
                  '  -u <num>     update visualization every <num> steps\n' +
                  '  -dt <num>    use <num> as length of timestep\n',
                  '  -e           normalize total energy in each timestep',
                  file=sys.stderr)
            return
        i += 2
    simulation(num_particles, steps, num_threads, normalize_energy, update_interval)
