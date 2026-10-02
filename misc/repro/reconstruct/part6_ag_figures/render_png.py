# Comparison tool for the item ag_figures (not part of the authors' package): renders the old and the regenerated
# fig/ag_almost.pdf and fig/ag_cancel.pdf to PNG with PyMuPDF (the poppler module of Hopper has no pdftoppm), whole pages
# and each panel stacked old above new, plus close views of the dips, for the comparison by eye.
# Usage: python3 render_png.py (old_ag_*.pdf and ag_*.pdf in the working directory); PyMuPDF 1.28.2 was installed with
# pip install --no-deps --target into the folder below, outside the shared venv.
import sys
sys.path.insert(0, '/scratch/vsokolov/rh_book_repro/agents/part6_ag_figures/pylib')
import pymupdf
from PIL import Image, ImageDraw


def render(fn, clip=None, dpi=150):
    page = pymupdf.open(fn)[0]
    pix = page.get_pixmap(dpi=dpi, clip=pymupdf.Rect(*clip) if clip else None)
    return Image.frombytes('RGB', (pix.width, pix.height), pix.samples)


def stack(fig, clip, dpi, out):
    ims = [render('old_%s.pdf' % fig, clip, dpi), render('%s.pdf' % fig, clip, dpi)]
    w = max(i.width for i in ims); h = sum(i.height for i in ims) + 60
    canvas = Image.new('RGB', (w, h), 'white'); d = ImageDraw.Draw(canvas); y = 0
    for lab, im in zip(('old', 'new'), ims):
        d.text((5, y + 5), lab, fill='black'); canvas.paste(im, (0, y + 20)); y += im.height + 30
    canvas.save(out)
    print(out, canvas.size)


for fig in ('ag_almost', 'ag_cancel'):
    for ver in ('old_', ''):
        im = render('%s%s.pdf' % (ver, fig), dpi=150); im.save('%s%s.png' % (ver, fig)); print('%s%s.png' % (ver, fig), im.size)
    for k, (x0, x1) in enumerate([(0, 312), (312, 604), (604, 900)]):   # the three panels, top to bottom of the page
        stack(fig, (x0, 0, x1, 252), 220, 'cmp_%s_panel%d.png' % (fig, k))
# close views of the dips (PDF points, origin at the top left): ag_almost at theta = 292.4, ag_cancel at theta = gamma_0
stack('ag_almost', (435, 0, 520, 252), 500, 'cmp_ag_almost_dip1.png')
stack('ag_almost', (725, 0, 810, 252), 500, 'cmp_ag_almost_dip2.png')
stack('ag_cancel', (428, 0, 510, 252), 500, 'cmp_ag_cancel_dip1.png')
stack('ag_cancel', (719, 0, 801, 252), 500, 'cmp_ag_cancel_dip2.png')
