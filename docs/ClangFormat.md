# Clang Format

`clang-format` is a set of tools to format your codebase based on given styles(or styles created by the user in a config file). It can be used for the following languages: C/C++/Java/JavaScript/JSON/Objective-C/Protobuf/C#. It is a project of LLVM. More in depth information about the actual documentation can be found [here](https://clang.llvm.org/docs/ClangFormat.html)

## Format using a bash script
In the current project, we have implemented two ways to format the code. The first way is a "manual" format. We created a bash script that will automatically format every `.cpp` and `.h` file presented inside the root directory of the repository. Unfortunately, `clang-format` does not have the recursive functionality and thus there is a need to find every file and act the formatting independently.

### 1. Finding the root directory

In order to find the root directory, and actually `cd` in it, we had to use the following command:

```
cd -- "$(find ~/ -type d -name ReCoDE-SPH-solver-2D-NS | head -1)"
```

This command tries to find a directory (`-type d`) with a name of the repository's name (`-name ReCoDE-SPH-solver-2D-NS`). Note that, in the case that the repository had a different name, this flag's value should have been changed as well. We then take this result of the `find` command and pipe(`|`) it to `head` which will give back the first(`-1`) line (this is done in case the given name is presented many times inside the machine). The result of this combined effort is then passed to `cd` so we can change the directory successfully.

### 2. Creating the style config

`clang-format` works either in default mode (with `llvm` as the style guide), with a given style (like the one we used, `google`) or with a user defined style config.

```
if ! [ -f "$config_file" ]; then
    clang-format --style=google --dump-config > "$config_file"
fi
```

This command should search the current directory for the specified config file and if it's not present create it on the spot to be used with `clang-format`.

### 3. Find and Format every CPP file

The final command should do the actual work. Now that we have find the proper directory we can move on to formatting.

```
find . \( -name '*.cpp' -o -name '*.h' \) -exec clang-format -i {} +
```

This command should find every `.cpp` and `.h` file inside the current directory and then pass them to `clang-format` to format it. If the find returns any CPP file then the formatting will be happening on these files, otherwise the command won't fail but the find will pass an empty string to `clang-format` (that's why you need `+`).

### Usage

In order to use the formatter you need to just type `./clang-format.sh`, given that the user has the right permissions to execute the file, and the formatter should do all the rest. If the user wants to change the style of the formatting, they need to navigate to `clang-format.sh` and change the `--style` flag's value to something different. Few already built suggestions can be found [here](https://clang.llvm.org/docs/ClangFormatStyleOptions.html).

## Format on commit

Apart from the "manual" format, we also utilized the `pre-commit` Python package to automate the process of formatting our staged files on every `git commit`. We have created a configuration file for `pre-commit` (`.pre-commit-config.yaml`), where a pre-commit hook is specified, that will be installed to our local repository the first time we run the command `pre-commit install`. Then, every time we run `git commit`, if there are any "bad formatted" staged files, the package will identify them, automatically format them based on the clang-format configuration we have specified, and the commit will fail with a message that prompts us to re-stage the changes that the package made and re-commit. In the case that all the staged files are clang-formatted, the commit will succeed with an appropriate message. More information on the `pre-commit package` can be found [here](https://pre-commit.com/).

## Format Checking Pipeline

Finally, we have also implemented a final, holistic check of our code's format. Using a Github Actions workflow (`format.yml`), we have developed a pipeline that is trigerred every time we push code changes or start a Pull Request (PR). The pipeline checks all the `.cpp` and `.h` files in the repository, using [this](https://github.com/jidicula/clang-format-action) clang-format Github Action. If all the files are formatted based on the clang-format configuration we have specified, the pipeline succeeds, else it fails. This pipeline, combined with the appopriate branch rules, guarantees that our Github repository is always well-formated.