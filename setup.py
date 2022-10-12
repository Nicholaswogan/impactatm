from skbuild import setup

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
    
setup(
    name="impactatm",
    packages=['impactatm'],
    python_requires='>=3.6',
    version="0.1.0",
    license="MIT",
    install_requires=['numpy','scipy','cantera'],
    author='Nicholas Wogan',
    author_email = 'nicholaswogan@gmail.com',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url = "https://github.com/Nicholaswogan/impactatm",
    include_package_data=True,
    cmake_args=['-DSKBUILD=ON','-CMAKE_BUILD_TYPE=Release']
)

