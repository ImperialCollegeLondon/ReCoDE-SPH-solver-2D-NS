import os
import pandas as pd
import matplotlib.pyplot as plt

# Set current directory as the starting location for operations in this file
script_path = os.path.abspath(__file__)
os.chdir(os.path.dirname(script_path))


OUTPUT_LOCATION = "../exec/output"
DOMAIN_FILE_LOCATION = "../exec/input/domain.txt"
FILES = ["initial-positions.csv", "final-positions.csv"]  # Path to target text file


def get_axes_limits():
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


if __name__ == "__main__":
    for file in FILES:
        # Read position values from .csv file
        position_data = pd.read_csv(f"{OUTPUT_LOCATION}/{file}")

        position_x = position_data["Position_X"]
        position_y = position_data["Position_Y"]

        x_min, x_max, y_min, y_max = get_axes_limits()
        plt.xlim(x_min, x_max)  # Set x-axis limits from x_min to x_max
        plt.ylim(y_min, y_max)  # Set y-axis limits from y_min to y_max

        plt.scatter(position_x, position_y)
        plt.title(f'{file.split("-")[0].capitalize()} Positions')

        plt.xlabel("Position x (m)")
        plt.ylabel("Position y (m)")

        plt.show()
