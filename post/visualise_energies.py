import os
import matplotlib.pyplot as plt

# Set current directory as the starting location for operations in this file
script_path = os.path.abspath(__file__)
os.chdir(os.path.dirname(script_path))

OUTPUT_LOCATION = "../exec/output"
FILENAME = "energies.txt"


def read_txt_file(file_path):
    # Read file
    with open(file_path, "r") as file:
        # Skip first line that holds the titles
        next(file)

        # Split line into columns
        data = [line.split() for line in file]

        col1, col2, col3, col4 = zip(*data)

        # Convert each column to float
        col1_values = [float(value) for value in col1]
        col2_values = [float(value) for value in col2]
        col3_values = [float(value) for value in col3]
        col4_values = [float(value) for value in col4]

    return col1_values, col2_values, col3_values, col4_values


if __name__ == "__main__":
    time, kinetic_energy, potential_energy, total_energy = read_txt_file(
        f"{OUTPUT_LOCATION}/{FILENAME}"
    )

    plt.plot(time, total_energy, label="Total Energy")
    plt.plot(time, kinetic_energy, label="Kinetic Energy")
    plt.plot(time, potential_energy, label="Potential Energy")

    plt.title("Energy")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")

    plt.legend()
    plt.show()
