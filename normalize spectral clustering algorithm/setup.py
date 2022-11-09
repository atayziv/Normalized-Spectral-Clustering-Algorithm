from setuptools import find_packages, setup, Extension
from pathlib import Path

SRC_PATH = Path("./")
sources = [str(file) for file in SRC_PATH.glob("*.c")]
module1 = Extension("spkm", sources=sources)

setup(
    name="spkm",
    version="1.0",
    ext_modules=[module1],
    install_requires=["invoke"],
    packages=find_packages(),
)
