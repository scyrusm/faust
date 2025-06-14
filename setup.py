#from distutils.core import setup
from setuptools import setup, find_packages
import os
import io

def _read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop('encoding', 'utf-8')
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text


def get_requirements(path):
    content = _read(path)
    return [
        req
        for req in content.split("\n")
        if req != "" and not (req.startswith("#") or req.startswith("-"))
    ]


install_requires = get_requirements('requirements.txt')
setup(
    name='faust',
    version='0.1dev',
    author='Samuel C. Markson',
    author_email='smarkson@alum.mit.edu',
    install_requires=install_requires, 
    extras_require={"GUI": ["PySide6", "matplotlib"]},
    packages=find_packages(),
    scripts=['scripts/faust'],
    license='BSD',
    description='for CRISPR screen analysis',
    include_package_data = True,
    long_description=open('README.rst').read(),
)

