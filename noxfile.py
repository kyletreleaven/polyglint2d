import nox

nox.options.default_venv_backend = "uv"


@nox.session
def test(session):
    try:
        toml = nox.project.load_toml("pyproject.toml")
        deps = toml["project"]["dependencies"]
        test_deps = toml["project"]["optional-dependencies"]["test"]

    except:
        import json
        print(json.dumps(toml, indent=2))
        raise

    session.install(*deps, *test_deps)  # but not the package itself...
    session.run("pytest", *(session.posargs or []))  # posargs for test filtering


@nox.session
def format(session):
    try:
        toml = nox.project.load_toml("pyproject.toml")
        package_deps = toml["project"]["optional-dependencies"]["package"]

    except:
        import json
        print(json.dumps(toml, indent=2))
        raise

    session.install(*package_deps)
    session.run("black", "src")  # TODO: Migrate options to config?


@nox.session
def dev(session):
    session.install("-e", ".")

    try:
        toml = nox.project.load_toml("pyproject.toml")
        dev_deps = toml["project"]["optional-dependencies"]["dev"]

    except:
        import json
        print(json.dumps(toml, indent=2))
        raise

    session.install(*dev_deps)
    session.run("ipython")


@nox.session
def run(session):
    session.install("-e", ".")

    try:
        toml = nox.project.load_toml("pyproject.toml")
        dev_deps = toml["project"]["optional-dependencies"]["dev"]

    except:
        import json
        print(json.dumps(toml, indent=2))
        raise

    session.install(*dev_deps)
    session.run("python", *(session.posargs or []))


@nox.session
def notebook(session):
    session.install("-e", ".")

    try:
        toml = nox.project.load_toml("pyproject.toml")
        dev_deps = toml["project"]["optional-dependencies"]["dev"]

    except:
        import json
        print(json.dumps(toml, indent=2))
        raise

    session.install(*dev_deps)
    session.run("jupyter", "notebook")


@nox.session
def docs(session):
    """

    Configuring Sphinx for document generation seems like a travesty.

    https://eikonomega.medium.com/getting-started-with-sphinx-autodoc-part-1-2cebbbca5365
    https://www.youtube.com/watch?v=KKfQnxQBoWE
    https://stackoverflow.com/questions/2701998/automatically-document-all-modules-recursively-with-sphinx-autodoc/62613202#62613202

    """
    session.install("-e", ".")

    try:
        toml = nox.project.load_toml("pyproject.toml")
        pkg_deps = toml["project"]["optional-dependencies"]["package"]

    except:
        import json
        print(json.dumps(toml, indent=2))
        raise

    session.install(*pkg_deps)
    session.run("python", "-m", "sphinx_autobuild", "docs", "html")


@nox.session
def nbenv(session):
    """

    # https://nbconvert.readthedocs.io/en/latest/usage.html#notebook-and-preprocessors

    e.g.,
    # Execute the notebook in-place, i.e., as implicit test step for doc prep!
    jupyter nbconvert --execute --inplace docs/basic-demo.ipynb

    # Then, for version control..

    # Note quite enough, leaves the output metadata.
    jupyter nbconvert --clear-output --inplace docs/basic-demo.ipynb

    # Clear output _and_ metadata!
    jupyter nbconvert --inplace \
        --ClearOutputPreprocessor.enabled=True --ClearMetadataPreprocessor.enabled=True \
        docs/basic-demo.ipynb

    """
    session.install("-e", ".")
    toml = nox.project.load_toml("pyproject.toml")
    pkg_deps = toml["project"]["optional-dependencies"]["package"]
    session.install(*pkg_deps)
    session.run(*session.posargs)
    # session.run("make", "-C", "docs", "html")
