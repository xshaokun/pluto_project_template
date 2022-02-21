# Project Name

Github is a great platform to do version control (make your modification tracable), idea management (issue/project), collaborations (fork & pull request) and sharing.

For accomplishing a vision "*Every part in a research work should be repeatable*", this project template tries to offer a strategy on how to manage those research projects of simulaton based on Pluto code.

The organization of this git template is consistant with the following **Pluto Project Management Specification**.

## Directories

* `code/`: Source code of simulations.
* `runs/`: Setup and parameters of different models
* `plotting/`: scripts for data analysis
* `assets/`: some essential data or something else

## Pluto Project Management Specification (PPMS)

### Introduction

* Sometimes it is hard, after a period of time, to remember what modifications were made and which file the modifications were based on.

* Sometimes you redo a simulation but find the results are different and don't know why.

* Sometimes some tests could mess up your directory and you don't want to remove them.

* Sometimes too many folders for parameter investigation make you crazy.

* Sometimes ...

So let's make some rules to manage a simulation project better !

### Version Control

The files that are added to version control should be at least **necessary to be tracked** or **universal** (see `.gitignore` for reference).

For example, object files produced by compling are not necessary to be tracked and not universal on other machine either.

Other case for example, the executable `pluto` is not universal, but it is necessary to be tracked for re-running the simulation for parameter investigation in case source code were modified. While it is preferred to ignore the executable in the directory of source code and copy it to the output directory.

### PLUTO Source Code

The directory of Pluto's source code should be named as `$PLUTO_DIR` as recommended by the official reference manual, locating at somewhere outside here.

Any modifications or implementations in souce code for the project should be stored in `code/`, which is also the directory where compiling was done.

#### Modifications on source code

Any modifications on the source code, including `definitions.h` and `init.c`, should be made in a new feature branch. Only after a period of developlemt is finished, this branch can be merged into `main` branch and an new release with an elaborate description is required.

The model or the parameter utilized in source code should be with a comment of its reference.

All test branches should be preserved and never be merged into `main` branch, therefore no need to release them.

### The directory of outputs

The directory of outputs, which can be specified by keyword `output_dir` in `pluto.ini`, should be seperated from the directory of source code.

The output directory must also contain the following files:

* *a symbolic link* `pluto.ini` to the corresponding parameter file `*.ini` if it is outside the output directory
* *a copy of executable* `pluto` in case a restart is required
* *the log files* which recording the information about the simulation, e.g. `dbl.out`, `pluto.0.log`.

A batch of simulations require a commit for recording the version of source code they used.

If there are multiple simulations, the output directories and parameter files of multiple simulations should be named differently and with a `README.md` in `runs/` for specification. (For reference, see the organization of `runs/` and the workflow of a shell script of [Pluto job submission](https://github.com/xshaokun/gadgets/blob/579ed26d41f3eff12c39529baf7c80b7429b1aac/bin/psub).)

The output datasets should not be tracked.

### Scripts of Analysis

All the scripts of plotting figures of the paper should be stored in a separated directory with a `README.md` for specification. Here is `plotting/`.

The images should not be tracked.

### Supporting data

All the supporting data for the project should be stored in a separated directory with a  `README.md` for specification. Here is `assets/`.
