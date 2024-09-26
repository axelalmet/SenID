from setuptools import setup, find_packages

setup(
    name="SenD",
    version="0.1.0",
    author="Axel Almet",
    author_email="axelalmet@gmail.com",
    description="A package to infer senescence-leaning subclusters and SASP-driven cell-cell communication",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/axelalmet/SenD",
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    install_requires=[
        'numpy',
        'scipy',
        'anndata',
        'pandas',
        'scanpy'
    ],
)