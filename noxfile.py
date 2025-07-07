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
def build(session):
    toml = nox.project.load_toml("pyproject.toml")
    build_deps = toml["project"]["optional-dependencies"]["build"]
    session.install(*build_deps)  # but not the package itself...
    session.run("python", "-m", "build_util", *(session.posargs or []))  # posargs for test filtering


@nox.session
def format(session):
    try:
        toml = nox.project.load_toml("pyproject.toml")
        package_deps = toml["project"]["optional-dependencies"]["format"]

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
    session.install("-e", ".[docs]")
    session.run("mkdocs", "build")


@nox.session
def nbenv(session):
    """

    # https://nbconvert.readthedocs.io/en/latest/usage.html#notebook-and-preprocessors

    e.g.,
    # Execute the notebook in-place, i.e., as implicit test step for doc prep!
    jupyter nbconvert --execute --inplace docs/basic-demo.ipynb

    # Then, for version control..

    # Note quite enough, leaves the output metadata.
    nox -r --session nbenv -- \
        jupyter nbconvert --clear-output --inplace docs/visualization.ipynb

    nox -r --session nbenv -- \
        jupyter nbconvert --inplace \
        --ClearOutputPreprocessor.enabled=True \
        demo/visualization.ipynb

    --ClearMetadataPreprocessor.enabled=True \

    """
    session.install("-e", ".")
    toml = nox.project.load_toml("pyproject.toml")
    pkg_deps = toml["project"]["optional-dependencies"]["docs"]
    session.install(*pkg_deps)
    session.run(*session.posargs)
    # session.run("make", "-C", "docs", "html")
