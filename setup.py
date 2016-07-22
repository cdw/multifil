from setuptools import setup

setup(name='multifil',
      version='0.2',
      description='A spatial half-sarcomere model and the means to run it',
      url='https://github.com/cdw/multifil',
      author='C David Williams',
      author_email='cdave@uw.edu',
      license='MIT',
      packages=['multifil'],
      install_requires=['numpy', 'boto']
     )
