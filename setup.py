from setuptools import setup

setup(
    name='ggoutlier',
    version='0.0.1',
    url='https://github.com/pktrigg/ggoutlier',
    author=(
        "pktrigg;"
    ),
    author_email=(
        "pktrigg@gmail.com;"
    ),
    description=(
        'Tool for detecting isolated or anomalous depth measurements in bathymetry data'
    ),
    entry_points={
        "gui_scripts": [],
        "console_scripts": [
            'ggoutlier = ggoutlier:main',
        ],
    },
    packages=[
        '.',
    ],
    zip_safe=False,
    package_data={
        "ggoutlier": [
            "GGOutlierGIS.png",
            "Guardian.png",
        ]
    },
    install_requires=[
        'reportlab',
        'pyproj',
        'pyshp',
        'scikit-learn',
        'rasterio',
        'numpy',
        'matplotlib'
    ],
)
