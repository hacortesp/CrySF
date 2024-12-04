from setuptools import setup, find_packages
import os

src_dir = os.path.abspath(os.path.dirname(__file__))

# Setting up
setup(
    name='CrySF',
    version='1.0.0',
    author='Henry Andres Cortes Paez',
    author_email='<hacortes@bcamath.org>',
    description= 'Data mining clustering of ionic trajectories to determine the crystallographic sites occupied by diffusive ions within crystalline materials',
    maintainer='Henry Andres Cortes Paez',
    maintainer_email="hacortes@bcamath.org",
    packages=find_packages(),
    install_requires=['numpy', 'MDanalysis', 'matplotlib', 'scikit-learn'],
    license='BCAM License',
    long_description=open('README.md').read(),
    python_requires='>=3.8',
)
