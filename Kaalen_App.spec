# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

a = Analysis(
    ['App_PP.py'],
    pathex=['.'],
    binaries=[],
    datas=[
        # remove resources_rc.py if it's a Python import
    ],
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

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    name='Kaalen_App',
    debug=False,
    strip=False,
    upx=False,    # disable UPX on GitHub runner unless you install it
    console=False,
    icon="icon.ico",   # recommended to ensure your EXE has an icon
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    name='Kaalen_App'
)
