# Python Particle Simulation

A particle physics simulation written in Python. Particles interact via contact forces and are visualized using Tkinter.

## Dependencies

All dependencies are Python standard library modules. The only system-level requirement is **Tkinter** (`python3-tk`), used for optional visualization. No third-party pip packages are needed.

## Docker

### Prerequisites

- [Docker](https://docs.docker.com/get-docker/) installed and running on your machine.

### Build the image

From the root of the repository (the directory containing the `Dockerfile`), run:

```bash
docker build -t particle-sim .
```

This will:
1. Pull the `python:3-slim` base image.
2. Install the `python3-tk` and `tk-dev` system packages required by the graphics module.
3. Copy all project files into the container under `/app`.

### Run the container

Run the simulation with the default settings (20 particles, 10,000 steps, no visualization):

```bash
docker run --rm particle-sim
```

#### Pass command-line options

The simulation accepts several flags that can be passed after the image name:

| Flag | Description |
|------|-------------|
| `-n <num>` | Number of particles to simulate (default: 20) |
| `-s <num>` | Number of timesteps to run (default: 10000) |
| `-dt <num>` | Length of each timestep (default: 0.00005) |
| `-v` / `-g` | Enable visualization (requires a display) |
| `-u <num>` | Update visualization every `<num>` steps |
| `-e` | Normalize total energy each timestep |

Example — simulate 50 particles for 5,000 steps:

```bash
docker run --rm particle-sim python verlet.py -n 50 -s 5000
```

### Run the unit tests

```bash
docker run --rm particle-sim python runtests.py
```