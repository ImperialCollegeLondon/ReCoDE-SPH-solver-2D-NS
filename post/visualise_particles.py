import os
import matplotlib.pyplot as plt

# Set current directory as the starting location for operations in this file
script_path = os.path.abspath(__file__)
os.chdir(os.path.dirname(script_path))

output_location = "../exec/output"
files = ["initial-positions.txt", "final-positions.txt"]


def read_txt_file(file_path):
    # Read file
    with open(file_path, "r") as file:
        # Skip first line that holds the titles
        next(file)

        # Split line into columns
        data = [line.split() for line in file]

        col1, col2 = zip(*data)

        # Convert each column to float
        col1_values = [float(value) for value in col1]
        col2_values = [float(value) for value in col2]

    return col1_values, col2_values


if __name__ == "__main__":
    for file in files:
        position_x, position_y = read_txt_file(f"{output_location}/{file}")

        plt.scatter(position_x, position_y)
        plt.title(f"{file.split('-')[0].capitalize()} Positions")

        plt.xlabel("Position x (m)")
        plt.ylabel("Position y (m)")

        plt.show()
