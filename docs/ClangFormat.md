# Clang Format

`clang-format` is a set of tools to format your codebase based on given styles(or styles created by the user in a config file). It can be used for the following languages: C/C++/Java/JavaScript/JSON/Objective-C/Protobuf/C#. It is a project of LLVM. More in depth information about the actual documentation can be found [here](https://clang.llvm.org/docs/ClangFormat.html).

In this project, we have implemented two ways to format C++ code - one manual and one using git hooks.

## Format using a bash script

This approach involves a "manual" formatting. We've developed a bash script designed to automatically format each .cpp and .h file located in the root directory of the repository. Since `clang-format` lacks inherent recursive capabilities, it becomes necessary to identify each file and apply formatting individually.

### 1. Find the root directory

To find the root directory and navigate (i.e. `cd`) to it, we use the following command:

```
cd -- "$(find ~/ -type d -name ReCoDE-SPH-solver-2D-NS | head -1)"
```

This command looks for a directory (`-type d`) with a name matching that of the repository (`-name ReCoDE-SPH-solver-2D-NS`). It is important to note that, should the repository have a different name, the value of this flag should be adjusted accordingly. Subsequently, we capture the output of the `find` command and pipe(`|`) it to `head`, which extracts the first(`-1`) line (a precaution in case the specified name appears multiple times on the host machine). The outcome of this combined effort is then passed to `cd` to navigate to the correct location.

### 2. Create the style config

`clang-format` works either in default mode (with `llvm` as the style guide), with a given style (like the one we use, i.e. `google`) or with a user defined style config.

```
if ! [ -f "$config_file" ]; then
    clang-format --style=google --dump-config > "$config_file"
fi
```

This command searches the current directory for the specified config file and, in the case it's not present, creates it on the spot to be used by `clang-format`.

### 3. Format CPP and header files

The final command performs the code formatting for every CPP and header file.

```
find . \( -name '*.cpp' -o -name '*.h' \) -exec clang-format -i {} +
```

This command finds every `.cpp` and `.h` file inside the current directory and passes them to `clang-format` to format them. Should the `find` command locate any CPP or header files, formatting will be applied to these files. In the absence of matching files, the command won't encounter an error, but it's essential to include the + symbol in the `find` command. This ensures that even if no files are found, an empty string is passed to clang-format as expected.

### Usage

In order to use the formatter, you need to run the following command in the terminal:

```
./clang-format.sh
```

The formatter then proceeds to format all the detected files using the style rules provided in the clang configuration file. If the user wants to change the style of the formatting applied, they need to edit `clang-format.sh` to change the `--style` flag's value to the desired style guide for clang follow. Some suggestions can be found [here](https://clang.llvm.org/docs/ClangFormatStyleOptions.html).

## Format on commit

Apart from performing "manual" formatting using the "clang-format.sh" script, we utilised the `pre-commit` Python package to automate the process of formatting our staged files on every `git commit`. This way of code formatting is independent of the manual way and you don't need to use both to format your code or for them to work properly.

There is a configuration file for `pre-commit` (`.pre-commit-config.yaml`), in which a pre-commit git hook is specified. To install this hook to a local git repository the following command is required:

```
pre-commit install
```

Then, every time the `git commit` command is executed, any "bad formatted" **staged** files are detected and formatted by `pre-commit`, based on the clang-format configuration that has been specified. Moreover, the commit fails with a message that prompts the user to stage the files that were formatted by the package and make a new commit. In the case that all the staged files are formatted according to the specified style rules, the commit will succeed with an appropriate message.

More information on the `pre-commit package` can be found [here](https://pre-commit.com/).

## Format Checking Pipeline

Furthermore, we have implemented a final, holistic check of our code's format. Using a Github Actions workflow (defined in `format.yml`), we have developed a pipeline that is trigerred every time a git push happens or a Pull Request (PR) is opened.

The pipeline checks all the `.cpp` and `.h` files in the repository, using the clang-format Github Action defined [here](https://github.com/jidicula/clang-format-action). If all the files are formatted based on the clang-format configuration that has been specified, the pipeline succeeds, else it fails. This pipeline, combined with the appopriate branch rules, guarantees that our Github repository is always well-formated.