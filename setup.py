from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='cthreepo',
    version='0.1.1',
    author='Vamsi Kodali',
    author_email='vkkodali@gmail.com',
    description='A python tool to interconvert seq-ids in gff3, gtf, bed and other files.',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url='http://github.com/vkkodali/cthreepo',
    license='MIT',
    packages=['cthreepo'],
    zip_safe=False,
    entry_points = {
        'console_scripts': ['cthreepo=cthreepo:main'],
    },
    package_data={'cthreepo': ['mapfiles/*.map']},
    install_requires=['requests']
)
