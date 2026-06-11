# view_capture

**Capture a molecule's on-screen orientation from ASE's x3d viewer — and render it from that exact angle.**

When you make figures for an article you usually `view(atoms, viewer='x3d')`, rotate the model by hand
to a nice angle, and then *guess* the `rotation` string (`'-80x,10y,0z'`) for the POV-Ray render — over
and over. `view_capture` removes the guessing: rotate in the viewer, click **📷 Capture orientation**,
and you get the exact ASE `rotation` to drop into your render.

It is **one file** (`view_capture.py`). Everything below uses just that file.

---

## Install — pick any one

**1. Copy-paste.** Paste the contents of `view_capture.py` into a notebook cell. Done.

**2. Download + import** (matches how these notebooks already fetch `.traj`/`.cif` data):

```python
!wget -q https://raw.githubusercontent.com/MineBonik/view_capture/main/view_capture.py
from view_capture import view_capture, CAPTURED_VIEWS, set_capture
```

**3. pip install** (uses the bundled `pyproject.toml`):

```python
!pip install git+https://github.com/MineBonik/view_capture
from view_capture import view_capture, CAPTURED_VIEWS, set_capture
```

> The course notebooks (`1.3`, `1.4`) already include a **bootstrap cell** that tries `import
> view_capture` and, if missing, downloads it from `VIEW_CAPTURE_URL`. Set that URL once after you host
> the file (or just keep `view_capture.py` next to the notebooks).

---

## Use it

```python
view_capture(atoms, name='Pt111')      # same as view(atoms, viewer='x3d') + a button
```

1. Rotate the model with the mouse to the perspective you want.
2. Click **📷 Capture orientation**. You get a box showing the `rotation` string and matrix, the string
   is copied to your clipboard, and (in Colab) it is stored in `CAPTURED_VIEWS['Pt111']`.
3. Render from that angle:

```python
generic_projection_settings['rotation'] = CAPTURED_VIEWS['Pt111']['rotation']   # exact 3x3 matrix
write('Pt111.pov', Pt111, **generic_projection_settings, povray_settings=povray_settings)
```

**Not on Colab?** The auto-store is Colab-only, but the click still shows the string and copies it to
the clipboard. Paste it straight into your settings:

```python
generic_projection_settings['rotation'] = '-80.00x,10.00y,0.00z'
# or, to use the same dict workflow:
set_capture('Pt111', '-80.00x,10.00y,0.00z')
generic_projection_settings['rotation'] = CAPTURED_VIEWS['Pt111']['rotation']
```

Tips:
- Use a **distinct `name=`** per model so captures don't overwrite each other.
- Quick sanity check: click the button *without* rotating → it reads `0.00x,0.00y,0.00z`.

---

## What it stores

`CAPTURED_VIEWS[name]` is a dict with:

| key            | meaning |
|----------------|---------|
| `rotation`     | 3×3 NumPy rotation matrix — exact, gimbal-free; pass directly to `write(..., rotation=...)` |
| `rotation_str` | the readable `'..x,..y,..z'` string (your existing notebook style) |

`rotation` reproduces the **angle** you set. Framing is auto-fit by ASE; to **zoom**, add
`auto_bbox_size` to your render (smaller = closer, e.g. `0.8`; larger = more whitespace). For full
control of zoom *and* pan, pass an explicit `bbox=(xlo, ylo, xhi, yhi)` in Ångströms. Pan is not
captured.

---

## Requirements / portability

Needs only **`ase`, `numpy`, `IPython`** — already present in these notebooks, **no extra installs**.
The interactive viewer loads the **x3dom** JavaScript from `x3dom.org` at view time, exactly like ASE's
own `view(atoms, viewer='x3d')`. **So if the x3d viewer works for you, this file works.** The only
unsupported case is a browser with no internet — a limitation inherited from ASE's viewer, not added
here.

Works in Google Colab (auto-store) and any Jupyter / JupyterLab / VS Code notebook (clipboard + on-screen
string + `set_capture`).

## How it works (one line)

The browser reports the live camera as X3DOM's `viewMatrix()`; its transpose **is** ASE's rotation
matrix (ASE builds rotation from `right/up/-look`), which `ase.utils.irotate` turns into the `'..x,..y,..z'`
string. Verified end-to-end against known camera orientations.

## License

MIT — see the header in `view_capture.py`.
