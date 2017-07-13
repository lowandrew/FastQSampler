from setuptools import setup

setup(
	name = "FastQSampler",
	version = "0.1.-",
	scripts = ["sampler.py","sampler_wrapper.py"],
	install_requires=["biopython"],
	author="Andrew Low",
	author_email="andrew.low@inspection.gc.ca",
)