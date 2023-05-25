try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='AcidifyME',
      version='0.01',
      description='ME model acid stress',
      author='Bing Du',
      author_email='',
      url='https://github.com/bdu91/acidify-ME.git',
      packages=['AcidifyME'],
      )
