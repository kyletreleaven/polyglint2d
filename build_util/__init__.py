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

    output_dir = DOCS_DIR / "notebooks"
    output_file = output_dir / "demo.html"

    def create(self):

        noxrun("nbenv", [
            *"jupyter nbconvert --to html --execute".split(),
            # "--show-config-json",
            "--output", self.output_file,
            self.notebook_file
        ])
