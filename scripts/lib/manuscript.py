"""manuscript.py -- figure finishing for scripts/manuscript_figures.ipynb (2026-11-22).

The analysis libraries draw each figure with its own text and sizes.  This module takes such a
figure after it is drawn and lets the notebook cell change it without touching the library:

    with ms.capture():                      # the library's plt.show() is held back
        fig = cmp.fig_strain_solvent(...)
    ms.finish(fig, name='fig1_strain_solvent', text=TEXT, size=(7.0, 4.2), fonts=dict(label=11))

Every piece of text of a figure has a key (ms.text_items / ms.print_text):

    suptitle                     the figure title
    fig.text.<j>                 text placed on the figure itself (notes under the axes)
    ax<i>.title | .title_left | .title_right | .xlabel | .ylabel
    ax<i>.legend.title | ax<i>.legend.<j>      legend entries, in drawing order
    ax<i>.text.<j>               annotations and text boxes inside the axes
    (ax<i> = fig.axes[i]; a colorbar is an axes of its own, its label is that axes' ylabel or xlabel)

In a TEXT dict: a string replaces the text; '' removes it; None, or leaving the key out, keeps
what the library wrote (use that for text that carries a fitted number, so it follows the data).

Global settings: ms.set_style(...) once at the top of the notebook (STYLE below); each ms.finish
call can override them for one figure.
"""
from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path

import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------------------------
# global settings (ms.set_style changes them; None = leave what the library drew)
# ---------------------------------------------------------------------------------------------
STYLE = dict(
    font_family=None,        # e.g. 'CMU Serif', 'Times New Roman', 'Helvetica'; None keeps the notebooks' font
    mathtext_fontset=None,   # 'cm' | 'stix' | 'dejavusans' | 'dejavuserif' | 'custom'
    font_scale=1.0,          # multiplies every font size the library used (applied before the fixed sizes below)
    size_suptitle=None,      # fixed font sizes in points, per kind of text; None = library size x font_scale
    size_title=None,
    size_label=None,
    size_tick=None,
    size_legend=None,
    size_text=None,          # annotations, text boxes, notes on the figure
    lw_scale=1.0,            # multiplies every line width (curves, not the axes frame)
    ms_scale=1.0,            # multiplies every marker size
    spine_lw=None,           # axes frame line width
    tick_direction=None,     # 'in' | 'out' | 'inout'
    grid=None,               # True / False forces the grid on / off on every axes; None keeps it
    suptitle=True,           # False drops every figure title (journal figures carry it in the caption)
    panel_titles=True,       # False drops every axes title
    dpi=300,                 # of the saved raster files
    formats=('png', 'pdf'),  # one file per format
    transparent=False,
    out_dir=None,            # where ms.finish saves; None = do not save
)


def set_style(**kw):
    """Change the global settings (keys of STYLE).  Returns the settings in force."""
    bad = [k for k in kw if k not in STYLE]
    if bad:
        raise KeyError(f'unknown style key(s) {bad}; known: {list(STYLE)}')
    STYLE.update(kw)
    if STYLE['font_family']:
        plt.rcParams['font.family'] = STYLE['font_family']
        plt.rcParams['mathtext.rm'] = STYLE['font_family']
    if STYLE['mathtext_fontset']:
        plt.rcParams['mathtext.fontset'] = STYLE['mathtext_fontset']
    return dict(STYLE)


@contextmanager
def capture():
    """Hold back plt.show() while a library figure function runs, so the figure can still be edited
    (the inline backend closes a figure once it is shown).  The library's own saved file is unaffected."""
    old = plt.show
    plt.show = lambda *a, **k: None
    try:
        yield
    finally:
        plt.show = old


# ---------------------------------------------------------------------------------------------
# the text of a figure
# ---------------------------------------------------------------------------------------------
def text_items(fig):
    """[(key, artist, kind)] for every piece of text of the figure; kind = suptitle | title | label | legend | text."""
    items = []
    st = getattr(fig, '_suptitle', None)
    if st is not None:
        items.append(('suptitle', st, 'suptitle'))
    for j, t in enumerate(fig.texts):
        if t is not st:
            items.append((f'fig.text.{j}', t, 'text'))
    for i, ax in enumerate(fig.axes):
        p = f'ax{i}'
        for nm, art in (('title', ax.title), ('title_left', ax._left_title), ('title_right', ax._right_title)):
            if art.get_text():
                items.append((f'{p}.{nm}', art, 'title'))
        for nm, art in (('xlabel', ax.xaxis.label), ('ylabel', ax.yaxis.label)):
            if art.get_text():
                items.append((f'{p}.{nm}', art, 'label'))
        leg = ax.get_legend()
        if leg is not None:
            if leg.get_title().get_text():
                items.append((f'{p}.legend.title', leg.get_title(), 'legend'))
            for j, t in enumerate(leg.get_texts()):
                items.append((f'{p}.legend.{j}', t, 'legend'))
        for j, t in enumerate(ax.texts):
            items.append((f'{p}.text.{j}', t, 'text'))
    return items


def text_dict(fig):
    """{key: current string} of the figure, in drawing order."""
    return {k: a.get_text() for k, a, _ in text_items(fig)}


def print_text(fig, name='TEXT'):
    """Print the figure's text as a dict literal to paste into a cell and edit."""
    print(format_text(fig, name))


def format_text(fig, name='TEXT'):
    rows = [f'{name} = {{']
    for k, v in text_dict(fig).items():
        rows.append(f'    {k!r}: {_lit(v)},')
    rows.append('}')
    return '\n'.join(rows)


def _one(s):
    """Literal of a string without newlines: raw when it has backslashes (LaTeX), so it reads as it is typed."""
    if '\\' in s and not s.endswith('\\'):
        q = '"' if '"' not in s else ("'" if "'" not in s else None)
        if q:
            return 'r' + q + s + q
    return repr(s)


def _lit(s):
    """A string literal that keeps backslashes readable (raw) and writes each line break as a separate '\\n' piece."""
    parts = s.split('\n')
    if len(parts) == 1:
        return _one(s)
    return '(' + " '\\n' ".join(_one(p) for p in parts) + ')'


# ---------------------------------------------------------------------------------------------
# finishing
# ---------------------------------------------------------------------------------------------
def finish(fig, name=None, text=None, size=None, fonts=None, limits=None, legends=None, style=None, save=None, show=True):
    """Apply the cell's edits and the global settings to a drawn figure, save it and show it.

    name     file stem under STYLE['out_dir'] (one file per STYLE['formats'])
    text     {key: string | '' | None} (see the module docstring; ms.print_text(fig) lists the keys)
    size     (width, height) in inches
    fonts    font sizes in points for this figure: dict(suptitle=, title=, label=, tick=, legend=, text=), or
             {key: size} for one piece of text (same keys as `text`)
    limits   {'ax0': dict(xlim=(a, b), ylim=(a, b), xscale='log', xticks=[...], yticks=[...])}
    legends  {'ax0': dict(loc='upper right', ncol=2, frameon=False, ...)}; {'ax0': False} removes that legend
    style    STYLE keys overridden for this figure only
    save     False = do not save this one; None = save when STYLE['out_dir'] is set and name is given
    """
    if fig is None:
        print('nothing to finish (the figure was skipped: see the message above)')
        return None
    S = dict(STYLE, **(style or {}))
    text, fonts = dict(text or {}), dict(fonts or {})
    kinds = ('suptitle', 'title', 'label', 'tick', 'legend', 'text')
    if size is not None:
        fig.set_size_inches(*size)

    # ---- legends first: a rebuilt legend has new text artists ----
    for key, opt in (legends or {}).items():
        ax = fig.axes[int(key[2:])]
        leg = ax.get_legend()
        if leg is None:
            continue
        if opt is False:
            leg.remove()
            continue
        handles = list(leg.legend_handles)
        labels = [t.get_text() for t in leg.get_texts()]
        kw = dict(fontsize=leg.get_texts()[0].get_fontsize() if leg.get_texts() else None, framealpha=leg.get_frame().get_alpha())
        kw.update(opt)
        leg.remove()
        ax.legend(handles, labels, **kw)

    items = text_items(fig)
    known = {k for k, _, _ in items}
    unknown = [k for k in list(text) + [k for k in fonts if k not in kinds] if k not in known]
    if unknown:
        print(f'note: no such text in this figure: {unknown}  (ms.print_text(fig) lists the keys)')
    for key, art, kind in items:
        if key in text and text[key] is not None:
            art.set_text(text[key])
        if kind == 'suptitle' and not S['suptitle']:
            art.set_text('')
        if kind == 'title' and not S['panel_titles'] and key not in text:
            art.set_text('')
        fs = fonts.get(key, fonts.get(kind, S.get(f'size_{kind}')))
        art.set_fontsize(fs if fs is not None else art.get_fontsize() * S['font_scale'])
        if S['font_family']:
            art.set_fontfamily(S['font_family'])

    for i, ax in enumerate(fig.axes):
        fs = fonts.get('tick', S['size_tick'])
        for lab in ax.get_xticklabels() + ax.get_yticklabels():
            lab.set_fontsize(fs if fs is not None else lab.get_fontsize() * S['font_scale'])
            if S['font_family']:
                lab.set_fontfamily(S['font_family'])
        for off in (ax.xaxis.get_offset_text(), ax.yaxis.get_offset_text()):
            off.set_fontsize(fs if fs is not None else off.get_fontsize() * S['font_scale'])
        if S['lw_scale'] != 1.0 or S['ms_scale'] != 1.0:
            for ln in ax.get_lines():
                ln.set_linewidth(ln.get_linewidth() * S['lw_scale'])
                ln.set_markersize(ln.get_markersize() * S['ms_scale'])
        if S['spine_lw'] is not None:
            for sp in ax.spines.values():
                sp.set_linewidth(S['spine_lw'])
        if S['tick_direction']:
            ax.tick_params(which='both', direction=S['tick_direction'])
        if S['grid'] is not None:
            ax.grid(S['grid'])
        opt = dict((limits or {}).get(f'ax{i}', {}))
        for k, fn in (('xscale', ax.set_xscale), ('yscale', ax.set_yscale), ('xlim', ax.set_xlim), ('ylim', ax.set_ylim),
                      ('xticks', ax.set_xticks), ('yticks', ax.set_yticks)):
            if k in opt:
                fn(opt.pop(k))
        if opt:
            ax.set(**opt)

    do_save = (S['out_dir'] is not None and name is not None) if save is None else bool(save)
    if do_save and name is not None and S['out_dir'] is not None:
        out = Path(S['out_dir'])
        out.mkdir(parents=True, exist_ok=True)
        for ext in S['formats']:
            fig.savefig(out / f'{name}.{ext}', dpi=S['dpi'], bbox_inches='tight', transparent=S['transparent'])
        print('saved', out / name, '(' + ', '.join(S['formats']) + ')')
    if show:
        try:
            from IPython.display import display
            display(fig)
        except Exception:
            pass
        plt.close(fig)
    return fig


# ---------------------------------------------------------------------------------------------
# figures drawn here (composites the deck libraries do not have)
# ---------------------------------------------------------------------------------------------
def fig_displacement_evolution(comp, perm=None, lvl=None, perm_fit='drive', figsize=None):
    """Displacement evolution and its poroelastic fit, one panel per experiment:
      (a) the free re-swelling after the uniaxial compression hold (the level's unload): u_z/L since the plates
          stopped against the material coordinate, dots = profiles, lines = the cos(k pi zeta) fit -> D_c;
      (b) the permeation run: u_z/L_0 since the pressure ramp, dots = profiles, lines = the consolidation fit
          (support face pinned, feed face free) -> D_c; dotted = the fitted steady state.  perm_fit 'unload' draws
          the permeation run's own re-swelling instead (runs with a Phase 2c unload).
    comp / perm: decks from cmp.load_deck.  Panel (b) is left out when no permeation deck is given or its network
    was frozen (no displacement to fit)."""
    import numpy as np
    from matplotlib.colors import Normalize
    import uniaxial as uni

    cfg, levels = comp['cfg'], comp['levels']
    pick = [L for L in levels if str(L['lvl']) == str(lvl)] if lvl is not None else [L for L in levels if (L.get('unload') or {}).get('prof')]
    L = (pick or levels)[0]
    Pf = (L.get('unload') or {}).get('prof')
    F = None
    if perm is not None:
        P = perm['P']
        F = (P.get('unload') or {}).get('Dc') if perm_fit == 'unload' else P.get('Dc')
        if F is None:
            print(f"permeation panel left out: the run has no displacement fit ({'no unload' if perm_fit == 'unload' else 'frozen network or no profiles'})")
    if Pf is None and F is None:
        print('displacement-evolution figure skipped (no re-swelling profiles and no permeation fit)')
        return None
    n = int(Pf is not None) + int(F is not None)
    fig, axes = plt.subplots(1, n, figsize=figsize or (9.0 * n + 1.0, 6.8), squeeze=False, constrained_layout=True)
    axes = list(axes[0])
    k = 0
    if Pf is not None:
        ax = axes[k]
        norm, cmap = uni._reswelling_profiles_panel(cfg, ax, Pf, letter='(a)' if n > 1 else '', primary=False)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02).set_label('time step')
        k += 1
    if F is not None:
        ax = axes[k]
        zff = np.linspace(0.0, 1.0, 400)
        norm = Normalize(vmin=F['ts'][F['early'][0]], vmax=F['ts'][F['early'][-1]])
        cmap = plt.cm.viridis
        for i in F['early']:
            c = cmap(norm(F['ts'][i]))
            ax.plot(F['zl'][i], F['yl'][i], 'o', color=c, ms=2.5, alpha=0.4)
            ax.plot(zff, F['u_model'](zff, F['t_lj'][i]), '-', color=c, lw=1.6)
        ax.plot(zff, F['u_inf'](zff), 'k:', lw=1.8, label=r'steady state ($t\to\infty$, fitted modes)')
        ax.set(xlabel=r'$\zeta$ (material coordinate; 0 = support face, 1 = feed face)',
               ylabel=r'$(u_z-u_0)/L_0$' if F.get('rigid') == 'contact' else r'$u_z/L_0$')
        ax.set_title(('(b) ' if n > 1 else '') + rf"permeation, {'re-swelling' if perm_fit == 'unload' else 'consolidation'} profiles: "
                     rf"$D_c={uni.sig(F['Dc'])}$, $R^2={uni.sig(F['R2'])}$, {F['n_modes']} modes", fontsize=13)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=10, loc='best')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02).set_label('time step')
    fig.suptitle('Displacement evolution and the poroelastic fit: re-swelling after uniaxial compression'
                 + (' and permeation' if F is not None else ''), fontsize=13, fontweight='bold')
    return fig
