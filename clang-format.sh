#!/usr/bin/bash

cd -- "$(find ~/ -type d -name ReCoDE-SPH-solver-2D-NS | head -1)"
config_file=".clang-style"


if ! [ -f "$config_file" ]; then
    clang-format --style=google --dump-config > "$config_file"
fi


find . \( -name '*.cpp' -o -name '*.hpp' \) -exec clang-format -i {} +

rm "$config_file"

echo "Formatting complete."