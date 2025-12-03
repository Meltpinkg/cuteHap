# coding=utf-8

from setuptools import setup, find_packages, Extension

with open('README.md') as f:
    readme = f.read()

ext_modules = [
    Extension(
        "cuteHap.genotype_improve",
        sources=[
            "src/cuteHap/genotype_improve.c",
        ],
        extra_compile_args=['-std=c99', '-O3'],
    ),
    Extension(
        "cuteHap.prove_func",
        sources=[
            "src/cuteHap/prove_func.c",
        ],
        extra_compile_args=['-std=c99', '-O3'],
    ),
]

setup(
    name = "cutehap",
    version = "1.0.2",
    description = "Haplotype-resolved genomic structural variation detection with cuteHap",
    author = "Shuqi Cao",
    author_email = "sqcao@stu.hit.edu.cn",
    url = "https://github.com/Meltpinkg/cuteHap",
    license = "MIT",
    packages = find_packages("src"),
    package_dir = {"": "src"},
    include_package_data = True,
    package_data={
        'cuteHap': [
            '*.pyx',
            '*.c',
            '*.py',
        ],
    },
    data_files = [("", ["LICENSE"])],
    scripts=['src/cuteHap/cuteHap'],
    # long_description = LONG_DESCRIPTION,
    long_description = readme,
    long_description_content_type = 'text/markdown',
    zip_safe = False,
    install_requires = ['scipy', 'pysam', 'Biopython', 'Cigar', 'numpy<2.0.0', 'Cython'],
    ext_modules=ext_modules
)
