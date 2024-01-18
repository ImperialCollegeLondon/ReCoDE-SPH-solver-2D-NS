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

    time = energy_data["t"]
    kinetic_energy = energy_data["Ek"]
    potential_energy = energy_data["Ep"]
    total_energy = energy_data["Etotal"]

    plt.plot(time, total_energy, label="Total Energy")
    plt.plot(time, kinetic_energy, label="Kinetic Energy")
    plt.plot(time, potential_energy, label="Potential Energy")

    plt.title("Energy")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")

    plt.legend()
    plt.show()
