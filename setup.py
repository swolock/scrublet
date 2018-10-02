from distutils.core import setup

setup(
    name = "scrublet",
    packages = ['scrublet'],
    package_dir={'': 'src'},
    version = '0.2',
    description = 'Doublet prediction in single-cell RNA-sequencing data',
    author = 'Samuel L. Wolock',
    author_email = 'swolock@g.harvard.edu',
    url = 'https://github.com/swolock/scrublet',
    download_url = 'https://github.com/swolock/scrublet/tarball/0.1',
    install_requires=['cython', 'numpy', 'scipy', 'scikit-learn', 'scikit-image', 'matplotlib', 'annoy', 'numba', 'pandas', 'umap-learn'],
    )
