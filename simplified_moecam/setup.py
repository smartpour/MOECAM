from setuptools import setup

classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
]

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="moecam",
    version="1.0.0",
    description="Multi-Objective Evolutionary Comparison and Analysis Module",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/smartpour/moecam",
    author='MOECAM Team',
    author_email='contact@moecam.org',
    license='MIT',
    py_modules=['moecam'],
    package_dir={'': 'src'},
    install_requires=['cffi>=1.0.0', 'numpy>=1.19.0', 'matplotlib>=3.3.0'],
    setup_requires=['cffi>=1.0.0'],
    cffi_modules=['./src/build_moecam.py:ffibuilder'],
    include_package_data=True,
    package_data={'':['tests/test_*.py']},
    classifiers=classifiers,
    python_requires=">=3.7",
)
