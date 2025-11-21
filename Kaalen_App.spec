# -*- mode: python ; coding: utf-8 -*-
import os
import sys

block_cipher = None

a = Analysis(
    ['App_PP.py'],
    pathex=['.'],
    binaries=[],
    datas=[
        ('resources_rc.py', '.'),  # This is correct
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
        'scipy._lib.array_api_compat.numpy',
        'scipy._lib.array_api_compat.numpy.fft',
        'scipy.linalg.cython_blas', 
        'scipy.linalg.cython_lapack',
        'scipy.optimize.minpack',
        
        # CRITICAL: Add these for Qt resources
        'PyQt5.Qt',
        'PyQt5.QtWidgets',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

# Add Qt plugins explicitly
from PyInstaller.utils.hooks import collect_qt_plugins

# Collect Qt plugins (important for icons and styles)
binaries = []
binaries.extend(collect_qt_plugins('platforms'))
binaries.extend(collect_qt_plugins('styles'))
binaries.extend(collect_qt_plugins('imageformats'))

a.binaries.extend(binaries)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(pyz,
          a.scripts, 
          a.binaries,
          a.zipfiles,
          a.datas,
          name='Kaalen_App',
          debug=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=False, 
          disable_windowed_traceback=False,
          target_arch=None,
          codesign_identity=None,
          entitlements_file=None,
          # Remove the external icon reference since you use resources
          icon=None)  # No external icon file

coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='Kaalen_App')
