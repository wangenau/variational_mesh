from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='variational_mesh',
    version='0.0.1',
    description='Variational Mesh',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://gitlab.com/wangenau/variational_mesh',
    author='Wanja Schulze',
    author_email='wangenau@protonmail.com',
    license='APACHE2.0',
    packages=find_packages(exclude=('tests')),
    install_requires=['matplotlib', 'numpy', 'pyscf==1.7.1'],
    python_requires='>=3',
    include_package_data=True,
    zip_safe=False,
)
