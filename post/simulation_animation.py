import math
from datetime import datetime

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import pandas as pd


# ======= Global Variables =======
DOMAIN_FILE_LOCATION = "../exec/input/domain.txt"
OUTPUT_LOCATION = "../exec/output"

SCATTER_PARTICLES = []  # List to store scatter positions
LINES_ENERGY = []  # List to store energy lines
LINE_TIMESTEP = None  # List to store timestep lines

ENERGY_TYPES = ["Ek", "Ep", "Etotal"]
ENERGY_LABELS = ["Kinetic Energy", "Potential Energy", "Total Energy"]
ENERGY_COLOURS = ["blue", "red", "purple"]
POSITIONS_DF = None
ENERGIES_DF = None
TOTAL_FRAMES = None
FRAME_DELAYS = []


def get_axes_limits():
    """
    Get the axes limits for the positions plot from the domain file
    """
    with open(DOMAIN_FILE_LOCATION) as f:
        lines = f.readlines()
        for line in lines:
            if "left_wall" in line:
                x_min = float(line.split()[2])
            if "right_wall" in line:
                x_max = float(line.split()[2])
            if "bottom_wall" in line:
                y_min = float(line.split()[2])
            if "top_wall" in line:
                y_max = float(line.split()[2])

    return x_min, x_max, y_min, y_max


def calculate_num_particles():
    """
    Dynamically determine the number of particles used in the SPH simulation
    """
    timestamps = POSITIONS_DF["Timestamp"]
    num_particles = None
    prev_time = None
    curr_time = None
    for timestamp in timestamps:
        curr_time = timestamp

        if prev_time is not None and prev_time != curr_time:
            break
        else:
            num_particles = num_particles + 1 if num_particles is not None else 1
            prev_time = curr_time

    return num_particles


def create_figure():
    global SCATTER_PARTICLES, LINES_ENERGY, LINE_TIMESTEP
    # Set global Plot styles
    plt.rcParams["grid.alpha"] = 0.4
    plt.rcParams["axes.titlesize"] = 16  # Set the title font size
    plt.rcParams["axes.titlepad"] = 20
    plt.rcParams["xtick.major.pad"] = 10  # Padding for major x-axis ticks
    plt.rcParams["ytick.major.pad"] = 8  # Padding for major y-axis ticks
    plt.rcParams["axes.titleweight"] = "bold"  # Set the title font weight
    plt.rcParams["axes.labelsize"] = 12  # Font size for both x and y labels
    plt.rcParams["axes.labelpad"] = 12  # Font size for both x and y labels
    plt.rcParams["axes.spines.top"] = False

    # Set up the figure and axes
    _, axs = plt.subplots(
        3,
        1,
        figsize=(8, 15),
        gridspec_kw={"width_ratios": [1], "height_ratios": [1.5, 1, 1], "hspace": 0.5},
    )

    # Initialize empty scatter plots for positions of multiple particles
    for i in range(num_particles):
        scatter_particle = axs[0].scatter([], [], s=10, color="#0077BE")
        SCATTER_PARTICLES.append(scatter_particle)

    # Add multiple energy lines to the bottom subplot
    for i, energy_label in enumerate(ENERGY_LABELS):
        (line_energy,) = axs[1].plot(
            [], [], lw=2, label=energy_label, color=ENERGY_COLOURS[i]
        )
        LINES_ENERGY.append(line_energy)

    LINE_TIMESTEP = (axs[2].plot([], [], lw=2, color="#40C9A2"))[0]

    # Set up the positions axis attributes
    x_min, x_max, y_min, y_max = get_axes_limits()
    axs[0].set_xlim(x_min, x_max)
    axs[0].set_ylim(y_min, y_max + 0.05 * y_max)
    axs[0].set_xlabel("Position X")
    axs[0].set_ylabel("Position Y")
    axs[0].set_title("SPH Particle Positions")
    axs[0].set_aspect("equal")

    # Set up the energy axis attributes
    simulation_duration = ENERGIES_DF["Timestamp"].iloc[-1]
    axs[1].set_xlim(
        0, simulation_duration + simulation_duration * 0.05
    )  # Time progresses along the x-axis
    axs[1].set_ylim(0, ENERGIES_DF["Etotal"].max() + ENERGIES_DF["Etotal"].max() * 0.05)
    axs[1].set_xlabel("Time (s)")
    axs[1].set_ylabel("Energy (J)")
    axs[1].set_title("Energy Progression")
    axs[1].grid()  # Add gridlines to the bottom subplot
    axs[1].spines["right"].set_visible(False)  # Hide right spine
    axs[1].legend()  # Show legend
    axs[1].set_aspect("auto")

    axs[2].set_xlim(
        0, simulation_duration + simulation_duration * 0.05
    )  # Time progresses along the x-axis
    axs[2].set_ylim(
        0, ENERGIES_DF["Timestep"].max() + ENERGIES_DF["Timestep"].max() * 0.05
    )
    axs[2].set_xlabel("Time (s)")
    axs[2].set_ylabel("Timestep (dt)")
    axs[2].set_title("Adaptive Timestep")
    axs[2].grid()  # Add gridlines to the bottom subplot
    axs[2].spines["right"].set_visible(False)
    axs[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # plt.show()


def update(frame):
    """
    Update the plot for each provided frame
    """
    progress = round((frame / TOTAL_FRAMES * (frame * 100 / TOTAL_FRAMES)), 2)
    print(f"Progress: {str(progress)}%")

    # Update positions
    adjusted_frame = frame * num_particles
    for i, scatter_particle in enumerate(SCATTER_PARTICLES):
        particle_positions = POSITIONS_DF.iloc[
            adjusted_frame : adjusted_frame + num_particles
        ]
        x = particle_positions["Position_X"]
        y = particle_positions["Position_Y"]
        scatter_particle.set_offsets(list(zip(x, y)))

    # Update energy
    timestamps = ENERGIES_DF["Timestamp"].iloc[: frame + 1]
    for i, energy_type in enumerate(ENERGY_TYPES):
        energy = ENERGIES_DF.iloc[: frame + 1][
            energy_type
        ]  # Update energy up to current frame
        LINES_ENERGY[i].set_data(timestamps, energy)

    LINE_TIMESTEP.set_data(timestamps, ENERGIES_DF.iloc[: frame + 1]["Timestep"])

    return *SCATTER_PARTICLES, *LINES_ENERGY, LINE_TIMESTEP


def create_animation():
    """
    Creates the animation of the paricles' positions and energies
    """
    global TOTAL_FRAMES
    # Determine the number of frames to display in the animation
    TOTAL_FRAMES = len(ENERGIES_DF)

    desired_animation_duration = math.ceil(
        ENERGIES_DF["Timestamp"].iloc[-1] + 1
    )  # Animation duration in seconds
    frame_rate = 30  # Frame rate of the animation in frames per second (fps)
    frames_progress = TOTAL_FRAMES // (desired_animation_duration * frame_rate)
    frames_to_display = range(0, TOTAL_FRAMES, frames_progress)

    # Build the animation
    ani = animation.FuncAnimation(
        plt.gcf(), update, frames=frames_to_display, blit=True, interval=10000
    )

    # Save the animation
    animation_filename = (
        f'fluid_simulation_{datetime.now().strftime("%Y%m%d_%H%M%S")}.mp4'
    )
    ani.save(animation_filename, fps=frame_rate)

    print("Progress 100%")
    print("Done!")
    print("The animation was saved to " + animation_filename)


if __name__ == "__main__":
    # Read position and energy output data
    POSITIONS_DF = pd.read_csv(
        f"{OUTPUT_LOCATION}/simulation-positions.csv"
    )  # Load the position data
    ENERGIES_DF = pd.read_csv(f"{OUTPUT_LOCATION}/energies.csv")  # Load the energy data

    # Calculate the number of particles in the simulation
    num_particles = calculate_num_particles()

    # Create the figure
    create_figure()

    # Create and store the animation
    create_animation()
