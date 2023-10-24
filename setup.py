from setuptools import find_packages, setup
import neokant


VERSION = neokant.VERSION


# parses requirements from file
with open("requirements.txt") as f:
    required = f.read().splitlines()

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# Build the Python package
setup(
    name='neokant',
    version=VERSION,
    packages=find_packages(exclude=["legacy"]),
    entry_points={
        'console_scripts': [
            'neokant-annotator=neokant.command_line:annotate',
            'neokant-parser=neokant.command_line:parse_output',
            'neokant-database=neokant.command_line:build_database',
        ],
    },
    author="TRON - Translational Oncology at the University Medical Center of the Johannes Gutenberg University Mainz"
    "- Computational Medicine group",
    author_email='johannes.hausmann@tron-mainz.de',
    description='Screen large collections of healthy tissue data for expression of neoantigen candidates',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    requires=[],
    install_requires=required,
    classifiers=[
        'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3 :: Only',
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix"
      ],
    python_requires='>=3.10',
    license='MIT'
)
