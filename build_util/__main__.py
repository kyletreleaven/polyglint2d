import click

from . import *

@click.command
@click.argument("target")
def cli(target: str):
    target = globals()[target]
    target.create()

cli()
