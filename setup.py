
from setuptools import setup, find_packages

setup(
	name = "polyglint2d",
	description = "A small utility for integration of linear functions over interiors of polygons in 2d",
	author = "Kyle Treleaven",
	author_email = "ktreleav@gmail.com",
	version = "0.0.0",
	install_requires = ['numpy', 'scipy']
	packages = find_packages(),
	namespace_packages = [ 'setiptah', ],
)

