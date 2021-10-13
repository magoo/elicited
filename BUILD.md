**Package / pypi Build**
- Empty the `dist/` directory
- `python3 -m build`
- `python3 -m twine upload dist/*`

**Rebuild from Pypi**
- `pip install --upgrade elicited --force-reinstall` 