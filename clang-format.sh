#!/usr/bin/bash

# Find the directory searching from users home directory
cd -- "$(find ~/ -type d -name ReCoDE-SPH-solver-2D-NS | head -1)"
config_file=".clang-format"

# Creating a style file based on a given format (we chose Google's)
if ! [ -f "$config_file" ]; then
    clang-format --style=google --dump-config > "$config_file"
fi

# Search for every .cpp and .h file and format them
find . \( -name '*.cpp' -o -name '*.h' \) -exec clang-format -i {} +

rm "$config_file"

echo "Formatting complete."