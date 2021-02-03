from setuptools import find_packages, setup

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='var_mesh',
    version='1.0.2',
    description='Variational Mesh',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://gitlab.com/wangenau/variational_mesh',
    author='Wanja Schulze',
    author_email='wangenau@protonmail.com',
    license='APACHE2.0',
    packages=find_packages(exclude=('docs', 'examples', 'public')),
    install_requires=['matplotlib', 'numpy', 'pyscf>=1.7'],
    python_requires='>=2.7 ,>=3.4',
    include_package_data=True,
    zip_safe=False,
)
