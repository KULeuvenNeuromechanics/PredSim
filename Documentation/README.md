
Documentation
=============

The PredSim documentation is available on [the website](https://kuleuvenneuromechanics.github.io/PredSim/).

The [/source](./source) folder contains the source files for the website.
Source files are written in [reStructuredText](https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html) (preferred) or [markdown](https://www.markdownguide.org/cheat-sheet/).

## Updating the website

Updating the website consists of 3 steps.
1. Update the source files
2. Build the html files
3. Deploy to github pages

### 1. Update the source files
Open the files in [/source](./source) in a text editor or IDE and edit them.
Mind the different syntax between [.rst](https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html) and [.md](https://www.markdownguide.org/cheat-sheet/) files.

### 2. Build the html files

The html files for the website can be build by running the [build_website workflow](../.github/workflows/build_website.yml) through [github actions](https://github.com/KULeuvenNeuromechanics/PredSim/actions/workflows/build_website.yml). This creates an *artifact* github-pages.zip that can be downloaded.
github-pages.zip contains artifact.tar, which contains the website files. Open index.html with a browser to preview the website.

The html files for the website can also be build locally.
First, set up a python environment, e.g. with anaconda.
- `conda create -n website`
- `conda activate website`
- `python -m pip install --upgrade pip`
- `cd c/path/to/PredSim/Documentation`
- `pip install -r requirements.txt`

After it is set up, access it via
- `conda activate website`
- `cd c/path/to/PredSim/Documentation`

To build the website files in *build/html*, run `sphinx-build -M html source/ build/`.

To build an interactive preview that updates when you change the source files, run `sphinx-autobuild source/ build/`.





### 3. Deploy to github pages

Pushing a change to the master branch (i.e. merging a pull request) will run the [deploy workflow](../.github/workflows/deploy.yml).
This workflow builds and deploys the website based on the updated source files in the master branch.




