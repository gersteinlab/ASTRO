from setuptools import setup, find_packages

setup(
    name="ASTRO",
    version="0.2",                   
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'ASTRO = ASTRO.ASTRO_run:main',
            'filtmatbyrt = ASTRO.ASTRO_run:filtmatbyrt'
        ]
    },  
    install_requires=[],             
)

