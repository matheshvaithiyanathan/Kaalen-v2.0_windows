# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

# Standard analysis
a = Analysis(
    ['App_PP.py'],
    pathex=['.'],
    binaries=[],
    datas=[],
    hiddenimports=[
        'matplotlib.backends.backend_qt5agg',
        'scipy.special.orthogonal',
        'numpy.core._dtype_ctypes',
        'pandas._libs.tslibs.timedeltas',
        'PyQt5.QtNetwork',
        'PyQt5.QtPrintSupport',
        'PyQt5.QtCore',
        'PyQt5.QtGui',
        'PyQt5.QtWidgets',
        'PyQt5.Qt',
        'scipy._lib.array_api_compat.numpy',
        'scipy._lib.array_api_compat.numpy.fft',
        'scipy.linalg.cython_blas',
        'scipy.linalg.cython_lapack',
        'scipy.optimize.minpack',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    cipher=block_cipher,
)

# Add Qt plugins properly
from PyInstaller.utils.hooks import collect_qt_plugins
qt_plugins = collect_qt_plugins('all')   # safer

a.binaries.extend(qt_plugins)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    [],
    a.datas,
    name='Kaalen_App',
    debug=False,
    strip=False,
    upx=True,
    console=False,
    icon=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
)
