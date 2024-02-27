import os
import pandas as pd
import matplotlib.pyplot as plt

# Set current directory as the starting location for operations in this file
script_path = os.path.abspath(__file__)
os.chdir(os.path.dirname(script_path))


OUTPUT_LOCATION = "../exec/output"
FILENAME = "energies.csv"


if __name__ == "__main__":
    # Read energy values from .csv file
    energy_data = pd.read_csv(f"{OUTPUT_LOCATION}/{FILENAME}")

    time = energy_data["Timestamp"]
    timestep = energy_data["Timestep"]
    kinetic_energy = energy_data["Ek"]
    potential_energy = energy_data["Ep"]
    total_energy = energy_data["Etotal"]

    _, (nrg_ax, dt_ax) = plt.subplots(
        1,
        2,
        figsize=(15, 5),
        gridspec_kw={"width_ratios": [1, 1], "height_ratios": [1]},
    )

    nrg_ax.plot(time, total_energy, label="Total Energy")
    nrg_ax.plot(time, kinetic_energy, label="Kinetic Energy")
    nrg_ax.plot(time, potential_energy, label="Potential Energy")

    nrg_ax.set_title("Energy")
    nrg_ax.set_xlabel("Time (s)")
    nrg_ax.set_ylabel("Energy (J)")

    nrg_ax.legend()

    dt_ax.plot(time, timestep, label="Adaptive Timestep")

    dt_ax.set_title("Adaptive Timestep")
    dt_ax.set_xlabel("Time (s)")
    dt_ax.set_ylabel("Timestep (dt)")

    plt.show()
