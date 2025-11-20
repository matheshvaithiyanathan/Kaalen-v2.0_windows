# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

# 1. Define the primary executable analysis
a = Analysis(
    ['App_PP.py'],
    pathex=['.'],  # The current directory containing App_PP.py
    binaries=[],
    datas=[
        ('resources_rc.py', '.'), # Include the resource file
    ],
    hiddenimports=[
        'matplotlib.backends.backend_qt5agg',
        'scipy.special.orthogonal', # Often needed by SciPy
        'numpy.core._dtype_ctypes', # General PyInstaller helper for NumPy
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

# 2. Add dependencies required by PyQt5, pandas, numpy, scipy, lmfit, pyqtgraph
# PyInstaller normally handles most of these via hooks, but explicitly confirming 
# non-Python data files often required by packages:
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)

# 3. Collect necessary DLLs and binary files 
# This gathers the C/C++ libraries used by PyQt5, numpy, scipy, etc.
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
          console=False, # Use False for a GUI application without a command window
          disable_window_shadow=False,
          target_arch=None,
          codesign_identity=None,
          entitlements_file=None ,)
