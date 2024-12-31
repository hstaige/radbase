from setuptools import setup

if __name__ == '__main__':
    setup(
        name="radbase",
        install_requires=['matplotlib', 'numpy', 'pandas', 'lmfit', 'uncertainties'],
        packages=['radbase']
    )
