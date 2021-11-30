See: https://packaging.python.org/tutorials/packaging-projects/a
**Local install**

- `pip install .` (inside `elicited` dir)

**Package / pypi Build**
- Empty the `dist/` directory
- `python3 -m build`
- `python3 -m twine upload dist/*`

**Rebuild from Pypi**
- `pip install --upgrade elicited --force-reinstall` 

If `No module named build`, install it: 

- `python3 -m pip install --upgrade build`
- `python3 -m pip install --upgrade twine`