# Use the official Python 3 slim image as the base
FROM python:3-slim

# Install python3-tk (Tkinter) required by graphics.py for visualization.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3-tk \
        tk-dev \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory inside the container
WORKDIR /app

# Copy all project files into the container
COPY . .

# Default command: run the particle simulation (no visualization by default)
CMD ["python", "verlet.py"]
