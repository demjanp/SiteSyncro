# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

a = Analysis(['bin\\process.py'],
			pathex=['.venv\\Lib\\site-packages'],
			binaries=[],
			datas=[('src\\sitesyncro\\graphviz', 'sitesyncro\\graphviz'),
				('src\\sitesyncro\\pygraphviz', 'sitesyncro\\pygraphviz')],
			hiddenimports=['sitesyncro', 'matplotlib', 'matplotlib.backends.backend_pdf'],
			hookspath=[],
			runtime_hooks=[],
			excludes=[],
			win_no_prefer_redirects=False,
			win_private_assemblies=False,
			cipher=block_cipher,
			noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='sitesyncro',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True,
          icon=None,
          onefile=True)  # Set onefile to True