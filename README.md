# morphogenDiffusion-ABM
The program simulate an agent-based model of morphogen diffusion in extracellular space on a two-dimesional grid.

Developer - Amit Landge (amit.landge@uni-konstanz.de)

Usage-
morphogenDiffusionABM.py [-h] [-n NUM [NUM ...]]
                                [-gS GRIDSIZE [GRIDSIZE ...]] [-st STEPS]
                                [-nrw NARROWINTERFACE]
                                [-slr SAILOR [SAILOR ...]] [-lv LIVE]
                                [-uImg USEIMG] [-iName IMGNAME]
                                [-parScr PARAMETERSCREEN] [-parID PARAMETERID]

Available user options.

optional arguments:
  - h, --help            show this help message and exit
  - n NUM [NUM ...], --num NUM [NUM ...]
                        Number of agents. default=1
  - gS GRIDSIZE [GRIDSIZE ...], --gridSize GRIDSIZE [GRIDSIZE ...]
                        Size of the 2D grid. default=[100,100]
  - st STEPS, --steps STEPS
                        Number of simulation steps. default=500
  - nrw NARROWINTERFACE, --narrowInterface NARROWINTERFACE
                        Steps to reduce interface width. default=5
  - slr SAILOR [SAILOR ...], --sailor SAILOR [SAILOR ...]
                        Types of sailor (list). default=0
  - lv LIVE, --live LIVE
                        Whether to show plots during runtime or not.
                        default=No
  - uImg USEIMG, --useImg USEIMG
                        Whether to use given image to make the grid.
                        default=No
  - iName IMGNAME, --imgName IMGNAME
                        Name of the image file. default=s8192_4_scaled.tiff
  - parScr PARAMETERSCREEN, --parameterScreen PARAMETERSCREEN
                        Whether to perform a parameter screen. default=No
  - parID PARAMETERID, --parameterID PARAMETERID
                        Row number in the parameter file. default=0

Requires- Python3
  - alabaster=0.7.12=py37_0
  - anyio=3.5.0=py37haa95532_0
  - appdirs=1.4.4=pyhd3eb1b0_0
  - apptools=5.1.0=pyhd3eb1b0_0
  - argon2-cffi=21.3.0=pyhd3eb1b0_0
  - argon2-cffi-bindings=21.2.0=py37h2bbff1b_0
  - arrow=0.13.1=py37_0
  - astroid=2.6.6=py37haa95532_0
  - atomicwrites=1.4.0=py_0
  - attrs=21.4.0=pyhd3eb1b0_0
  - autopep8=1.6.0=pyhd3eb1b0_0
  - babel=2.9.1=pyhd3eb1b0_0
  - backcall=0.2.0=pyhd3eb1b0_0
  - bcrypt=3.2.0=py37he774522_0
  - binaryornot=0.4.4=pyhd3eb1b0_1
  - black=19.10b0=py_0
  - bleach=4.1.0=pyhd3eb1b0_0
  - brotli=1.0.9=ha925a31_2
  - brotlipy=0.7.0=py37h2bbff1b_1003
  - bzip2=1.0.8=he774522_0
  - ca-certificates=2022.2.1=haa95532_0
  - certifi=2021.10.8=py37haa95532_2
  - cffi=1.15.0=py37h2bbff1b_1
  - chardet=4.0.0=py37haa95532_1003
  - charset-normalizer=2.0.4=pyhd3eb1b0_0
  - click=8.0.4=py37haa95532_0
  - cloudpickle=2.0.0=pyhd3eb1b0_0
  - colorama=0.4.4=pyhd3eb1b0_0
  - configobj=5.0.6=py37haa95532_1
  - cookiecutter=1.7.2=pyhd3eb1b0_0
  - cryptography=36.0.0=py37h21b164f_0
  - curl=7.80.0=h2bbff1b_0
  - cycler=0.11.0=pyhd3eb1b0_0
  - debugpy=1.5.1=py37hd77b12b_0
  - decorator=5.1.1=pyhd3eb1b0_0
  - defusedxml=0.7.1=pyhd3eb1b0_0
  - diff-match-patch=20200713=pyhd3eb1b0_0
  - docutils=0.17.1=py37haa95532_1
  - entrypoints=0.3=py37_0
  - envisage=6.0.1=pyhd3eb1b0_0
  - expat=2.4.4=h6c2663c_0
  - fipy=3.4.2.1=py37h03978a9_2
  - flake8=3.9.2=pyhd3eb1b0_0
  - fonttools=4.25.0=pyhd3eb1b0_0
  - freetype=2.10.4=hd328e21_0
  - future=0.18.2=py37_1
  - hdf4=4.2.13=h712560f_2
  - hdf5=1.10.6=nompi_h5268f04_1114
  - icc_rt=2019.0.0=h0cc432a_1
  - icu=58.2=ha925a31_3
  - idna=3.3=pyhd3eb1b0_0
  - imagesize=1.3.0=pyhd3eb1b0_0
  - importlib-metadata=4.8.2=py37haa95532_0
  - importlib_metadata=4.8.2=hd3eb1b0_0
  - inflection=0.5.1=py37haa95532_0
  - intel-openmp=2022.0.0=haa95532_3663
  - intervaltree=3.1.0=pyhd3eb1b0_0
  - ipykernel=6.9.1=py37haa95532_0
  - ipython=7.31.1=py37haa95532_0
  - ipython_genutils=0.2.0=pyhd3eb1b0_1
  - isort=5.9.3=pyhd3eb1b0_0
  - jedi=0.18.1=py37haa95532_1
  - jinja2=2.11.3=pyhd3eb1b0_0
  - jinja2-time=0.2.0=pyhd3eb1b0_2
  - joblib=1.1.0=pyhd3eb1b0_0
  - jpeg=9d=h2bbff1b_0
  - json5=0.9.6=pyhd3eb1b0_0
  - jsoncpp=1.8.4=h1ad3211_1002
  - jsonschema=3.2.0=pyhd3eb1b0_2
  - jupyter_client=6.1.12=pyhd3eb1b0_0
  - jupyter_core=4.9.2=py37haa95532_0
  - jupyter_server=1.13.5=pyhd3eb1b0_0
  - jupyterlab=3.2.9=pyhd8ed1ab_0
  - jupyterlab_pygments=0.1.2=py_0
  - jupyterlab_server=2.10.3=pyhd3eb1b0_1
  - keyring=23.4.0=py37haa95532_0
  - kiwisolver=1.3.2=py37hd77b12b_0
  - lazy-object-proxy=1.6.0=py37h2bbff1b_0
  - libblas=3.9.0=8_mkl
  - libcblas=3.9.0=8_mkl
  - libcurl=7.80.0=h86230a5_0
  - libiconv=1.15=h1df5818_7
  - liblapack=3.9.0=8_mkl
  - libnetcdf=4.7.4=nompi_h3a9aa94_107
  - libpng=1.6.37=h2a8f88b_0
  - libspatialindex=1.9.3=h6c2663c_0
  - libssh2=1.9.0=h7a1dbc1_1
  - libtiff=4.2.0=hd0e1b90_0
  - libwebp=1.2.2=h2bbff1b_0
  - libxml2=2.9.12=h0ad7f3c_0
  - lz4-c=1.9.2=h62dcd97_2
  - m2w64-gcc-libgfortran=5.3.0=6
  - m2w64-gcc-libs=5.3.0=7
  - m2w64-gcc-libs-core=5.3.0=7
  - m2w64-gmp=6.1.0=2
  - m2w64-libwinpthread-git=5.0.0.4634.697f757=2
  - mako=1.1.4=pyhd3eb1b0_0
  - markdown=3.3.4=py37haa95532_0
  - markupsafe=1.1.1=py37hfa6e2cd_1
  - matplotlib-base=3.5.1=py37hd77b12b_0
  - matplotlib-inline=0.1.2=pyhd3eb1b0_2
  - mayavi=4.7.1=py37h11817fe_2
  - mccabe=0.6.1=py37_1
  - mistune=0.8.4=py37hfa6e2cd_1001
  - mkl=2020.4=hb70f87d_311
  - mpmath=1.2.1=py37haa95532_0
  - msys2-conda-epoch=20160418=1
  - munkres=1.1.4=py_0
  - mypy_extensions=0.4.3=py37haa95532_1
  - nbclassic=0.3.5=pyhd3eb1b0_0
  - nbclient=0.5.11=pyhd3eb1b0_0
  - nbconvert=6.1.0=py37haa95532_0
  - nbformat=5.1.3=pyhd3eb1b0_0
  - nest-asyncio=1.5.1=pyhd3eb1b0_0
  - notebook=6.4.8=py37haa95532_0
  - numpy=1.21.5=py37h5fa1a60_0
  - numpydoc=1.2=pyhd3eb1b0_0
  - openssl=1.1.1m=h2bbff1b_0
  - packaging=21.3=pyhd3eb1b0_0
  - pandas=1.3.5=py37h9386db6_0
  - pandocfilters=1.5.0=pyhd3eb1b0_0
  - paramiko=2.8.1=pyhd3eb1b0_0
  - parso=0.8.3=pyhd3eb1b0_0
  - pathspec=0.7.0=py_0
  - patsy=0.5.2=py37haa95532_1
  - pdoc3=0.9.2=pyhd3eb1b0_0
  - pexpect=4.8.0=pyhd3eb1b0_3
  - pickleshare=0.7.5=pyhd3eb1b0_1003
  - pillow=9.0.1=py37hdc2b20a_0
  - pip=21.2.4=py37haa95532_0
  - pluggy=1.0.0=py37haa95532_0
  - poyo=0.5.0=pyhd3eb1b0_0
  - prometheus_client=0.13.1=pyhd3eb1b0_0
  - prompt-toolkit=3.0.20=pyhd3eb1b0_0
  - psutil=5.8.0=py37h2bbff1b_1
  - ptyprocess=0.7.0=pyhd3eb1b0_2
  - pycodestyle=2.7.0=pyhd3eb1b0_0
  - pycparser=2.21=pyhd3eb1b0_0
  - pydocstyle=6.1.1=pyhd3eb1b0_0
  - pyface=7.3.0=py37haa95532_1
  - pyflakes=2.3.1=pyhd3eb1b0_0
  - pygments=2.11.2=pyhd3eb1b0_0
  - pylint=2.9.6=py37haa95532_1
  - pyls-spyder=0.4.0=pyhd3eb1b0_0
  - pynacl=1.4.0=py37h62dcd97_1
  - pyopenssl=22.0.0=pyhd3eb1b0_0
  - pyparsing=3.0.4=pyhd3eb1b0_0
  - pyqt=5.9.2=py37hd77b12b_6
  - pyrsistent=0.18.0=py37h196d8e1_0
  - pysocks=1.7.1=py37_1
  - python=3.7.12=h7840368_100_cpython
  - python-dateutil=2.8.2=pyhd3eb1b0_0
  - python-lsp-black=1.0.0=pyhd3eb1b0_0
  - python-lsp-jsonrpc=1.0.0=pyhd3eb1b0_0
  - python-lsp-server=1.2.4=pyhd3eb1b0_0
  - python-slugify=5.0.2=pyhd3eb1b0_0
  - python_abi=3.7=2_cp37m
  - pytz=2021.3=pyhd3eb1b0_0
  - pywin32=302=py37h827c3e9_1
  - pywin32-ctypes=0.2.0=py37_1001
  - pywinpty=2.0.2=py37h5da7b33_0
  - pyyaml=6.0=py37h2bbff1b_1
  - pyzmq=22.3.0=py37hd77b12b_2
  - qdarkstyle=3.0.2=pyhd3eb1b0_0
  - qstylizer=0.1.10=pyhd3eb1b0_0
  - qt=5.9.7=vc14h73c81de_0
  - qtawesome=1.0.3=pyhd3eb1b0_0
  - qtconsole=5.2.2=pyhd3eb1b0_0
  - qtpy=1.11.2=pyhd3eb1b0_0
  - regex=2021.8.3=py37h2bbff1b_0
  - requests=2.27.1=pyhd3eb1b0_0
  - rope=0.22.0=pyhd3eb1b0_0
  - rtree=0.9.7=py37h2eaa2aa_1
  - scikit-learn=1.0.2=py37hcabfae0_0
  - scipy=1.7.3=py37hb6553fb_0
  - seaborn=0.11.2=hd8ed1ab_0
  - seaborn-base=0.11.2=pyhd8ed1ab_0
  - send2trash=1.8.0=pyhd3eb1b0_1
  - setuptools=58.0.4=py37haa95532_0
  - sip=4.19.13=py37hd77b12b_0
  - six=1.16.0=pyhd3eb1b0_1
  - sniffio=1.2.0=py37haa95532_1
  - snowballstemmer=2.2.0=pyhd3eb1b0_0
  - sortedcontainers=2.4.0=pyhd3eb1b0_0
  - sphinx=4.4.0=pyh6c4a22f_1
  - sphinxcontrib-applehelp=1.0.2=pyhd3eb1b0_0
  - sphinxcontrib-devhelp=1.0.2=pyhd3eb1b0_0
  - sphinxcontrib-htmlhelp=2.0.0=pyhd3eb1b0_0
  - sphinxcontrib-jsmath=1.0.1=pyhd3eb1b0_0
  - sphinxcontrib-qthelp=1.0.3=pyhd3eb1b0_0
  - sphinxcontrib-serializinghtml=1.1.5=pyhd3eb1b0_0
  - spyder=5.1.5=py37haa95532_1
  - spyder-kernels=2.1.3=py37haa95532_0
  - sqlite=3.38.0=h2bbff1b_0
  - statsmodels=0.13.0=py37h2bbff1b_0
  - sympy=1.9=py37h03978a9_1
  - tbb=2020.2=h2d74725_4
  - terminado=0.13.1=py37haa95532_0
  - testpath=0.5.0=pyhd3eb1b0_0
  - text-unidecode=1.3=pyhd3eb1b0_0
  - textdistance=4.2.1=pyhd3eb1b0_0
  - threadpoolctl=2.2.0=pyh0d69192_0
  - three-merge=0.1.1=pyhd3eb1b0_0
  - tinycss=0.4=pyhd3eb1b0_1002
  - tk=8.6.11=h2bbff1b_0
  - toml=0.10.2=pyhd3eb1b0_0
  - tornado=6.1=py37h2bbff1b_0
  - traitlets=5.1.1=pyhd3eb1b0_0
  - traits=6.2.0=py37h2bbff1b_0
  - traitsui=7.2.1=pyhd3eb1b0_0
  - typed-ast=1.4.3=py37h2bbff1b_1
  - typing-extensions=3.10.0.2=hd3eb1b0_0
  - typing_extensions=3.10.0.2=pyh06a4308_0
  - ujson=4.0.2=py37hd77b12b_0
  - unidecode=1.2.0=pyhd3eb1b0_0
  - urllib3=1.26.8=pyhd3eb1b0_0
  - vc=14.2=h21ff451_1
  - vs2015_runtime=14.27.29016=h5e58377_2
  - vtk=8.2.0=py37h359a0f6_218
  - watchdog=2.1.6=py37haa95532_0
  - wcwidth=0.2.5=pyhd3eb1b0_0
  - webencodings=0.5.1=py37_1
  - websocket-client=0.58.0=py37haa95532_4
  - wheel=0.37.1=pyhd3eb1b0_0
  - whichcraft=0.6.1=pyhd3eb1b0_0
  - win_inet_pton=1.1.0=py37haa95532_0
  - wincertstore=0.2=py37haa95532_2
  - winpty=0.4.3=4
  - wrapt=1.12.1=py37he774522_1
  - xz=5.2.5=h62dcd97_0
  - yaml=0.2.5=he774522_0
  - yapf=0.31.0=pyhd3eb1b0_0
  - zipp=3.7.0=pyhd3eb1b0_0
  - zlib=1.2.11=h8cc25b3_4
  - zstd=1.4.8=h498f385_0
