import nox

nox.options.default_venv_backend = "uv"


@nox.session
def test(session):
    """Runs the unit tests."""
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
    """Minimal session to run build commands.

    e.g., `nox -r -s build -- {{target}}`

    """
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
    """Opens an IPython session with the package installed."""
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
    """A session that hosts the repo's notebooks for development."""
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
    """Compiles the mkdocs documentation.

    See `build_util` for dependencies.

    """
    session.install("-e", ".[docs]")
    session.run("mkdocs", "build")


@nox.session
def nbenv(session):
    """Use this session to run nbconvert utils.

    # https://nbconvert.readthedocs.io/en/latest/usage.html#notebook-and-preprocessors

    """
    session.install("-e", ".")
    toml = nox.project.load_toml("pyproject.toml")
    pkg_deps = toml["project"]["optional-dependencies"]["docs"]
    session.install(*pkg_deps)
    session.run(*session.posargs)
