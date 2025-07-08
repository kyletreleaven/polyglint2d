import subprocess as sp
import pathlib

DIR = pathlib.Path(__file__).parent
PROJECT_DIR = DIR.parent
BUILD_DIR = PROJECT_DIR / "build"
DOCS_DIR = PROJECT_DIR / "docs"

def singleton(factory):
    return factory()


def noxrun(session, args_, *args, **kwargs):
    return sp.run([
        "nox", "-r", "-s", session, "--", *args_
    ], *args, **kwargs)


@singleton
class docs:
    """Compiles documentation including a Jupyter notebook."""

    def create(self):

        demo_notebook.create()  # TODO: Skip if it already exists.

        noxrun("docs", [])


@singleton
class demo_notebook:
    """A nice notebook.

    https://nbconvert.readthedocs.io/en/latest/usage.html#

    """
    DIR = PROJECT_DIR / "demo"
    notebook_file =  DIR / "demo.ipynb"

    notebooks_dir = DOCS_DIR / "notebooks"
    images_dir = notebooks_dir / "images"
    output_file = notebooks_dir / "demo.html"

    def create(self):

        preprocs = ["build_util.extract_notebook_gif.ExtractGifPreprocessor"]
        incl_preproc = f"--Exporter.preprocessors={repr(preprocs)}"

        noxrun("nbenv", [
            *"jupyter nbconvert --to html --execute".split(),
            # "--show-config-json",
            incl_preproc,
            "--output-dir", self.notebooks_dir,
            "--output", self.output_file,
            self.notebook_file
        ])
