# coding=utf-8

from setuptools import setup, find_packages
from Cython.Build import cythonize

with open('README.md') as f:
    readme = f.read()

cython_extensions = cythonize("src/cuteHap/*.pyx")

setup(
    name = "cutehap",
    version = "1.0.1",
    description = "Haplotype-resolved genomic structural variation detection with cuteHap",
    author = "Shuqi Cao",
    author_email = "sqcao@stu.hit.edu.cn",
    url = "https://github.com/Meltpinkg/cuteHap",
    license = "MIT",
    packages = find_packages("src"),
    package_dir = {"": "src"},
    data_files = [("", ["LICENSE"])],
    scripts=['src/cuteHap/cuteHap'],
    # long_description = LONG_DESCRIPTION,
    long_description = readme,
    long_description_content_type = 'text/markdown',
    zip_safe = False,
    setup_requires = ['Cython'],
    install_requires = ['scipy', 'pysam', 'Biopython', 'Cigar', 'numpy', 'Cython'],
    ext_modules=cython_extensions
)
