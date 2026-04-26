from setuptools import setup, find_packages

setup(
    name="ASTRO",
    version="1.0",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "ASTRO = ASTRO.ASTRO_run:main",
            "ASTROutils = ASTRO.ASTRO_utils:main",
            "filtmatbyrt = ASTRO.ASTRO_run:filtmatbyrt",
        ]
    },
    extras_require={
        "velocity": [
            "anndata>=0.8",
            "matplotlib>=3.5",
            "numpy<2",
            "pandas>=1.3",
            "scanpy>=1.9",
            "scipy>=1.7",
            "scvelo>=0.2.5",
        ],
    },
    install_requires=[],
)
