from setuptools import setup

setup(name='cthreepo',
    version='0.1',
    description='A python tool to interconvert seq-ids in gff3, gtf, bed and other files.',
    url='http://github.com/vkkodali/cthreepo',
    author='Vamsi Kodali',
    author_email='vkkodali@gmail.com',
    license='MIT',
    packages=['cthreepo'],
    zip_safe=False,
    entry_points = {
        'console_scripts': ['cthreepo=cthreepo:main'],
    },
    package_data={'cthreepo': ['mapfiles/*.map']},
    install_requires=['requests']
)
