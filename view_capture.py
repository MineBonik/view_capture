# view_capture.py
#
# MIT License
# Copyright (c) 2026 ViLab Tartu (electrochemistry group) and contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
"""Capture a molecule's on-screen orientation from ASE's x3d viewer.

`view_capture(atoms, name=...)` is a drop-in replacement for
`view(atoms, viewer='x3d')` that adds a **Capture orientation** button. Rotate
the model in the viewer to the perspective you want for a figure, click the
button, and the live camera angle is converted to ASE's `rotation` --- both a
3x3 matrix and the `'..x,..y,..z'` string used by POV-Ray / matplotlib renders.
No more guessing Euler angles by trial and error.

How the captured value reaches Python
-------------------------------------
* In **Google Colab** the click stores it in ``CAPTURED_VIEWS[name]`` automatically.
* **Everywhere** the rotation string is shown on screen and copied to the clipboard.
* In any environment you can also type it back in with
  ``set_capture(name, '-80x,10y,0z')`` (paste the string the button gave you).

Then render with the captured angle::

    generic_projection_settings['rotation'] = CAPTURED_VIEWS['Pt111']['rotation']
    write('Pt111.pov', Pt111, **generic_projection_settings, povray_settings=povray_settings)

To zoom the rendered image, add ``auto_bbox_size`` to the render (smaller = closer,
e.g. ``0.8``); ASE auto-fits the framing otherwise.

Three ways to use this one file
-------------------------------
1. **Copy-paste** its contents into a notebook cell.
2. **Download + import** (matches how the course notebooks fetch data files)::

       !wget -q https://raw.githubusercontent.com/MineBonik/view_capture/main/view_capture.py
       from view_capture import view_capture, CAPTURED_VIEWS, set_capture

3. **pip install** (with the accompanying ``pyproject.toml``)::

       !pip install git+https://github.com/MineBonik/view_capture
       from view_capture import view_capture, CAPTURED_VIEWS, set_capture

Requirements
------------
Only ``ase``, ``numpy`` and ``IPython`` --- all already present in these
notebooks. The interactive viewer loads the x3dom JavaScript library from
``x3dom.org`` at view time, exactly like ASE's own ``view(atoms, viewer='x3d')``.
So if the x3d viewer works for you, this file works too. (Truly offline use ---
no browser internet --- is the one unsupported case, inherited from ASE.)
"""
from string import Template
from uuid import uuid4

import numpy as np
from ase.io.x3d import pretty_print, x3d_atoms

__version__ = "0.1.0"
__all__ = ["view_capture", "build_html", "set_capture", "CAPTURED_VIEWS", "__version__"]

# Filled in by the JS -> kernel callback (Colab) or by set_capture(). Maps
# name -> dict with keys: 'rotation' (3x3 np.ndarray) and 'rotation_str' (str).
CAPTURED_VIEWS = {}


def _fmt(a):
    """Format an angle to 2 dp, turning -0.00 into 0.00."""
    a = round(float(a), 2)
    if a == 0:
        a = 0.0
    return f'{a:.2f}'


def _rotation_str(R):
    """ASE 'Ax,Ay,Az' string reconstructing the rotation matrix R (via irotate)."""
    from ase.utils import irotate
    x, y, z = irotate(np.asarray(R, dtype=float))
    return f'{_fmt(x)}x,{_fmt(y)}y,{_fmt(z)}z'


def set_capture(name, rotation):
    """Manually store an orientation in CAPTURED_VIEWS (no viewer needed).

    Use this when you copied the rotation string from the button (or already
    know an angle) and want the same dict-based workflow as the auto-capture::

        set_capture('Pt111', '-80x,10y,0z')
        generic_projection_settings['rotation'] = CAPTURED_VIEWS['Pt111']['rotation']

    Args:
        name:     key to store under in CAPTURED_VIEWS.
        rotation: an ASE rotation string like '-80x,10y,0z', or a 3x3 matrix.

    Returns the stored entry.
    """
    from ase.utils import rotate
    if isinstance(rotation, str):
        R = rotate(rotation)
    else:
        R = np.asarray(rotation, dtype=float).reshape(3, 3)
    CAPTURED_VIEWS[name] = {
        'rotation': R,
        'rotation_str': _rotation_str(R),
    }
    return CAPTURED_VIEWS[name]


def _register_colab_callback():
    """Register the kernel callback once (no-op outside Colab)."""
    try:
        from google.colab import output as _colab_output
    except Exception:
        return False

    def _capture(name, mat9):
        R = np.array(mat9, dtype=float).reshape(3, 3)
        CAPTURED_VIEWS[name] = {
            'rotation': R,
            'rotation_str': _rotation_str(R),
        }
        return {'rotation_str': CAPTURED_VIEWS[name]['rotation_str']}

    _colab_output.register_callback('viewcapture.capture', _capture)
    return True


_HTML = Template(r"""
<link rel="stylesheet" type="text/css" href="https://www.x3dom.org/release/x3dom.css"></link>
<script type="text/javascript" src="https://www.x3dom.org/release/x3dom.js"></script>

<X3D id="x3d_$uid" showStat="false" showLog="false"
     style="width:$width; height:$height; border:1px solid #ccc;">
$scene
</X3D>

<div style="margin-top:6px; font-family:monospace;">
  <button id="btn_$uid" type="button" onclick="capture_$uid()"
          style="font-size:14px; padding:4px 10px; cursor:pointer;">
    &#128247; Capture orientation
  </button>
  <span id="msg_$uid" style="margin-left:8px; color:#2a7;"></span><br>
  <textarea id="out_$uid" rows="4" cols="78" readonly
            style="margin-top:6px; font-family:monospace; font-size:12px;"></textarea>
  <textarea id="clip_$uid" style="position:absolute; left:-9999px;"></textarea>
</div>

<script type="text/javascript">
(function() {
  function el(V, i, j) {
    return (typeof V.at === 'function') ? V.at(i, j) : V['_' + i + '' + j];
  }
  function sgn(x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); }
  function givens(a, b) {
    var c, s, r, cot, tan, u;
    if (b === 0) { c = sgn(a); s = 0; r = Math.abs(a); }
    else if (Math.abs(b) >= Math.abs(a)) {
      cot = a / b; u = sgn(b) * Math.sqrt(1 + cot * cot);
      s = 1 / u; c = s * cot; r = b * u;
    } else {
      tan = b / a; u = sgn(a) * Math.sqrt(1 + tan * tan);
      c = 1 / u; s = c * tan; r = a * u;
    }
    return [c, s, r];
  }
  function irotate(a) {
    var gx = givens(a[2][2], a[1][2]); var cx = gx[0], sx = gx[1], rx = gx[2];
    var gy = givens(rx, a[0][2]);      var cy = gy[0], sy = gy[1];
    var gz = givens(cx * a[1][1] - sx * a[2][1],
                    cy * a[0][1] - sy * (sx * a[1][1] + cx * a[2][1]));
    var cz = gz[0], sz = gz[1];
    return [Math.atan2(sx, cx) * 180 / Math.PI,
            Math.atan2(-sy, cy) * 180 / Math.PI,
            Math.atan2(sz, cz) * 180 / Math.PI];
  }
  function round6(v) { return Math.round(v * 1e6) / 1e6; }

  window.capture_$uid = function() {
    var x3dEl = document.getElementById('x3d_$uid');
    var out = document.getElementById('out_$uid');
    var msg = document.getElementById('msg_$uid');
    var rt = x3dEl ? x3dEl.runtime : null;
    if (!rt || typeof rt.viewMatrix !== 'function') {
      out.value = 'x3dom is not ready yet - rotate the model once, then click again.';
      return;
    }
    var V = rt.viewMatrix();  // world -> eye
    // R_ase columns = (right, up, -look) = rows of the view-matrix rotation,
    // i.e. R_ase = transpose(upper-left 3x3 of V).
    var R = [[el(V, 0, 0), el(V, 1, 0), el(V, 2, 0)],
             [el(V, 0, 1), el(V, 1, 1), el(V, 2, 1)],
             [el(V, 0, 2), el(V, 1, 2), el(V, 2, 2)]];
    var ang = irotate(R);
    function fa(v) { v = Math.round(v * 100) / 100; if (v === 0) v = 0; return v.toFixed(2); }
    var rs = fa(ang[0]) + 'x,' + fa(ang[1]) + 'y,' + fa(ang[2]) + 'z';

    var Rstr = '[[' + round6(R[0][0]) + ', ' + round6(R[0][1]) + ', ' + round6(R[0][2]) + '], ['
                    + round6(R[1][0]) + ', ' + round6(R[1][1]) + ', ' + round6(R[1][2]) + '], ['
                    + round6(R[2][0]) + ', ' + round6(R[2][1]) + ', ' + round6(R[2][2]) + ']]';
    out.value = "rotation string : " + rs + "\n"
              + "rotation matrix : " + Rstr + "\n"
              + "name            : '$name'\n"
              + "to zoom the render: add  auto_bbox_size=0.8  (smaller = closer)";

    try {
      var clip = document.getElementById('clip_$uid');
      clip.value = rs; clip.select(); document.execCommand('copy');
      msg.style.color = '#2a7';
      msg.textContent = '✓ rotation string copied to clipboard';
    } catch (e) {
      msg.style.color = '#a72';
      msg.textContent = '(copy blocked - select the text manually)';
    }

    if (window.google && google.colab && google.colab.kernel) {
      try {
        var flat = [R[0][0], R[0][1], R[0][2],
                    R[1][0], R[1][1], R[1][2],
                    R[2][0], R[2][1], R[2][2]];
        google.colab.kernel.invokeFunction('viewcapture.capture', ['$name', flat], {})
          .then(function(res) {
            try {
              var s = res.data['application/json'].rotation_str;
              msg.textContent = '✓ stored in CAPTURED_VIEWS["$name"]  (rotation=' + s + ')';
            } catch (e) {}
          });
      } catch (e) {}
    }
  };

  // x3dom auto-inits on page load; if output was injected late, reload once.
  window.addEventListener('load', function() {
    var elx = document.getElementById('x3d_$uid');
    if (window.x3dom && x3dom.reload && elx && !elx.runtime) {
      try { x3dom.reload(); } catch (e) {}
    }
  });
})();
</script>
""")


def build_html(atoms, name='view', width='400px', height='300px', uid=None):
    """Return the raw HTML string for the capturing viewer (no IPython needed)."""
    scene = pretty_print(x3d_atoms(atoms))
    return _HTML.safe_substitute(
        uid=uid or uuid4().hex[:8],
        name=str(name),
        width=width,
        height=height,
        scene=scene,
    )


def view_capture(atoms, name='view', width='400px', height='300px'):
    """Like `view(atoms, viewer='x3d')` but with a "Capture orientation" button.

    Args:
        atoms:  ase.Atoms to display.
        name:   key under which the captured orientation is stored in
                CAPTURED_VIEWS (Colab) and shown on screen. Use a distinct name
                per model so captures don't overwrite each other.
        width, height: CSS size of the viewer.

    Returns an IPython HTML object (display it as the last line of a cell).
    """
    from IPython.display import HTML

    _register_colab_callback()
    return HTML(build_html(atoms, name=name, width=width, height=height))
