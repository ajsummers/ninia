import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ninia",
    version="0.0.2",
    author="Alex Summers",
    author_email="ajs0201@auburn.edu",
    description="A small Python wrapper for Quantum Espresso - still in development",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ajsummers/ninia",
    project_urls={
        'Issues': 'https://github.com/ajsummers/ninia/issues'
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    include_package_data=True,
    package_data={
        'ninia': ['*.i', '*.sh']
    },
    python_requires=">=3.6",
)
