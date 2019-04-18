import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "scrublet",
    packages = ['scrublet'],
    package_dir={'': 'src'},
    version = '0.2',
    description = 'Doublet prediction in single-cell RNA-sequencing data',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author = 'Samuel L. Wolock',
    author_email = 'swolock@g.harvard.edu',
    url = 'https://github.com/allonkleinlab/scrublet',
    install_requires=['cython', 'numpy', 'scipy', 'scikit-learn', 'scikit-image', 'matplotlib', 'annoy', 'numba', 'pandas', 'umap-learn'],
    )
