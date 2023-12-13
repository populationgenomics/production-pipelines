# Requirements

You will need the following software installed to complete this setup:

- Homebrew (if you are on macOS), or an equivalent package manager
- Python 3.10.x (there have been some issues with 3.11)
- Java 8 or Java 11 (required for running Hail locally)

Since most of our org use Apple laptops for development, this rest of this guide assumes
that you are operating in a macOS environment.

## Homebrew

Set up homebrew according to the instructions [here](https://brew.sh/). Before closing
your terminal, please read any output from the installation script because it may print
additional instructions on how to complete the setup. Before proceeding, ensure that it
has installed correctly by running `brew doctor` in your terminal and checking that the
output does not contain any issues.

## Python

It is recommended that you install Python with pyenv, or some flavour of conda. This
document will guide you through installing alternative Python versions using a Python
version manager called [pyenv](https://github.com/pyenv/pyenv). You can skip this section
if you are using [conda](https://docs.conda.io/en/latest/) to manage your Python versions
and virtual environments. However, please make sure to create a new virtual environment
with **Python 3.10**, for example:

```sh
conda create --name pipelines python=3.10.12
```

### Pyenv Users

Take some time to understand how pyenv works by reading the short
[How It Works](https://github.com/pyenv/pyenv#how-it-works) section; this understanding
will help you to debug operating system and path related issues (i.e _'command not
found'_, or seemingly missing packages) when they occur. When you're ready, please
follow the setup instructions [here](https://github.com/pyenv/pyenv#installation).

Once you're done, open a new terminal and verify your pyenv installation by entering
`pyenv --version`. Now install Python 3.10 using

```sh
pyenv install 3.10.12
```

You may use a different 3.10 version which you can find by executing the command
`pyenv install --list` and scrolling to the top of the output. The rest of this guide
assumes that you are using `3.10.12`, so change any relevant commands to suit your
version of Python.

## Java

If using homebrew (recommended), install Java 11 using the command:

```sh
brew install openjdk@11
```

Follow the instructions that homebrew gives at the end of the installation which will
put Java into the system path so that other programs can find it. At the time of writing
this document, the command to do this is:

```sh
sudo ln -sfn \
    $HOMEBREW_PREFIX/opt/openjdk@11/libexec/openjdk.jdk \
    /Library/Java/JavaVirtualMachines/openjdk-11.jdk
```

**IMPORTANT:** this command might have changed, so it's important to follow the
instructions that homebrew gives you when running the install command above. Open a new
terminal and check that everything has installed correctly by running `java --version`.
You should see 11 as the major version in the output.

# Repository Setup

Open a new terminal and change directory to where your git repositories are stored. Now,
clone the repository recursively and `cd` into the folder:

```sh
git clone --recurse-submodules git@github.com:populationgenomics/production-pipelines.git
cd production-pipelines
```

### Conda Users

If using `conda`, activate the virtual environment that you created earlier.

### Pyenv Users

If using `pyenv`, execute `pyenv local 3.10.12` in your terminal. This will create an
application-specific file called `.python-version` that tells `pyenv` which version of
Python to use when inside this directory. Verify that you are using the correct version
of Python by running `python --version` in your terminal. Now create and activate a virtual
environment called `venv`:

```sh
python -m venv venv && source venv/bin/activate
```

Verify that you are using your virtual environment by confirming that the command
`which python` is pointing to the Python binary in the local `venv` folder that
was just created.

## Installation

Run the following pip commands:

```sh
pip install gnomad_methods/ .
pip install -r seqr-loading-pipelines/luigi_pipeline/requirements.txt
pip install .
pip install -r requirements-dev.txt
```

# Visual Studio Code

## Selecting a Python Interpreter

Press `cmd + shift + p` to open the command palette and type in `Select Interpreter`.
Click the menu item that says `Python: Select Interpreter`. Select the interpreter
option which has the path `./venv/bin/python`. VSCode might highlight this option by
placing a star to the left of it, and the text `(Recommended)` to the right. If using
conda, select the menu item relating to your virtual environment created earlier. See
Microsoft's documentation for more [details](https://code.visualstudio.com/docs/python/environments#_creating-environments).

## Setting up Testing

Set up the test runner in Visual Studio Code by following Microsoft's
[instructions](https://code.visualstudio.com/docs/python/testing#_configure-tests). When
asked, select `pytest` as your testing framework and `test` as your testing directory.
VSCode should now auto-detect tests and populate the sidebar with all the tests. Press
the `play` button to run the tests and make sure they're all passing. To debug tests,
place a breakpoint somewhere in the code, right-click on any test and select the
`Debug Test` menu item. Please use the interactive debugger to debug the application as
opposed to print statements as it is much more time efficient; being able to inspect
variables and step through code in real time is extremely helpful in understanding
complex applications.
