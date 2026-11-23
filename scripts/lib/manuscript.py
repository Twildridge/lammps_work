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

Figures drawn here rather than in a deck library: fig_displacement_evolution, and fig_stress_profiles
(solvent, polymer and total stress profiles across the cell, with the reservoirs, plates and pistons named),
and fig_gel_volume (with load_gel_volume: the Rg and BB gel-volume curves of an equilibration run).
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



# ---------------------------------------------------------------------------------------------
# partial-stress profiles across the cell (two-piston decks)
# ---------------------------------------------------------------------------------------------
SPECIES_COLORS = dict(s='#0072B2', p='#D55E00', t='#222222')   # solvent / polymer / total: Okabe-Ito blue and vermilion + black,
                                                               # distinct under every colour-vision deficiency and in grey
RESERVOIR_SHADE = dict(color='#56B4E9', alpha=0.10, lw=0, zorder=0)


def _fade(base, f):
    """Colour `base` mixed towards white: f = 0 white, 1 the colour itself."""
    from matplotlib.colors import to_rgb
    r, g, b = to_rgb(base)
    return (1 - f * (1 - r), 1 - f * (1 - g), 1 - f * (1 - b))


def fig_stress_profiles(deck, lvl=None, phases=('hold',), comp='zz', n_curves=None, figsize=None):
    """Solvent, polymer and total stress profiles across the cell and their evolution, with the regions named:
    feed reservoir | gel | permeate reservoir (each reservoir shaded from its piston to the gel face), and the feed
    piston, load plate, support plate and permeate piston drawn at their positions at the end of the phase.

    deck      a two-piston deck from cmp.load_deck (compression: one panel per entry of `phases`)
    phases    any of 'hold' (reference -> compressed, plates fixed) and 'unload' (compressed -> re-swollen)
    comp      'zz' | 'xx' | 'yy'
    n_curves  profiles drawn per species over the phase, evenly spaced in time (None = the deck Config's n_curves);
              early ones pale and thin, the last one bold.  Dashed = the eps = 0 reference.
    x is 1 - z/L_z as in the deck's own figures, feed side on the left; limits are just outside the two pistons."""
    import numpy as np
    from matplotlib.lines import Line2D
    import uniaxial as uni

    cfg, R = deck['cfg'], deck['R']
    L = deck['P'] if deck.get('kind') == 'permeation' else None
    if L is None:
        levels = deck['levels']
        L = ([x for x in levels if str(x['lvl']) == str(lvl)] or levels)[0]
    if comp not in L['stress'] or comp not in R['stress']:
        print(f'stress-profile figure skipped (no sigma_{comp} profiles loaded)')
        return None
    panels = []
    for ph in phases:
        D = L.get('U') if ph == 'unload' else L
        if D is None:
            print(f"panel '{ph}' left out: this level has no unload")
            continue
        panels.append((ph, D))
    if not panels:
        return None
    if n_curves is not None:
        cfg = type('C', (), dict(n_curves=int(n_curves)))()        # post_halt reads only n_curves
    n = len(panels)
    fig, axes = plt.subplots(n, 1, figsize=figsize or (13.0, 6.6 * n + 0.6), squeeze=False, constrained_layout=True)
    zx = uni.zn(R, R['z'])
    names = (('s', 'solvent', 's,'), ('p', 'polymer', 'p,'), ('t', 'total', ''))
    for k, (ph, D) in enumerate(panels):
        ax = axes[k, 0]
        S, Rs = D['stress'][comp], R['stress'][comp]
        t_end = float(np.asarray(D['ts'])[-1])
        # ---- the four boundaries at the end of the phase, in plot coordinates ----
        wz = D.get('wetz') or {}
        z_feed = wz['feed_at'](t_end) if callable(wz.get('feed_at')) else R.get('z_feed', np.nan)
        z_perm = wz['perm_at'](t_end) if callable(wz.get('perm_at')) else R.get('z_perm', np.nan)
        xb = {nm: float(uni.zn(R, z)) for nm, z in (('feed', z_feed), ('load', D.get('z_pist', R['z_piston'])),
                                                     ('supp', D.get('z_supp', R['z_support'])), ('perm', z_perm))}
        x_gel = sorted(float(uni.zn(R, z)) for z in (D['z_mem_lo'], D['z_mem_hi']))
        left_is_feed = xb['feed'] < xb['perm']
        ax.axvspan(*x_gel, **uni.GEL_SHADE)
        g_feed, g_perm = (x_gel[0], x_gel[1]) if left_is_feed else (x_gel[1], x_gel[0])     # the gel face on each side
        ax.axvspan(*sorted((xb['feed'], g_feed)), **RESERVOIR_SHADE)                        # solvent only: piston to gel face
        ax.axvspan(*sorted((g_perm, xb['perm'])), **RESERVOIR_SHADE)                        # (the plates stand in it once retracted)
        for nm, ls in (('feed', ':'), ('perm', ':'), ('load', '-.'), ('supp', '-')):
            if np.isfinite(xb[nm]):
                ax.axvline(xb[nm], color='0.15', ls=ls, lw=1.6, zorder=4)
        # ---- curves ----
        top = 0.0
        for key, _, _ in names:
            base = SPECIES_COLORS[key]
            ax.plot(zx, Rs[key + '_m'], '--', color=base, lw=1.8, zorder=3)
            ts, ev = uni.post_halt(cfg, D, D['ts'], S[key])
            m = len(ts)
            for i in range(m):
                last = i == m - 1
                ax.plot(zx, ev[i], '-', color=base if last else _fade(base, 0.25 + 0.45 * i / max(m - 1, 1)),
                        lw=3.0 if last else 1.2, zorder=6 if last else 5)
            inside = (zx > min(x_gel)) & (zx < max(x_gel))
            top = max(top, float(np.nanmax(ev[-1][inside])), float(np.nanmax(np.asarray(Rs[key + '_m'])[inside])))
        ax.axhline(0, color='k', lw=0.8, alpha=0.5, zorder=1)
        pad = 0.02
        lo, hi = min(xb['feed'], xb['perm']) - pad, max(xb['feed'], xb['perm']) + pad
        ax.set_xlim(lo, hi)
        ax.set_ylim(-0.04 * top, 1.30 * top)
        # ---- region and boundary names (text order: feed reservoir, gel, permeate reservoir, then the four boundaries) ----
        tr = ax.get_xaxis_transform()
        for x, s in ((0.5 * (xb['feed'] + g_feed), 'feed\nreservoir'), (0.5 * sum(x_gel), 'gel'),
                     (0.5 * (g_perm + xb['perm']), 'permeate\nreservoir')):
            ax.text(x, 0.965, s, transform=tr, ha='center', va='top', fontsize=19, fontweight='bold', linespacing=1.0, zorder=8)
        for nm, s in (('feed', 'feed piston'), ('load', 'load plate'), ('supp', 'support plate'), ('perm', 'permeate piston')):
            right = (xb[nm] > 0.5 * (lo + hi)) if nm in ('feed', 'perm') else (nm == 'supp') == left_is_feed      # plate names on the reservoir side
            ax.annotate(s, (xb[nm], 0.30), xycoords=tr, xytext=(7 if right else -7, 0), textcoords='offset points', rotation=90,
                        ha='left' if right else 'right', va='center', fontsize=15, color='0.15', zorder=8)
        ax.set_xlabel(r'position across the cell,  $1-z/L_z$' if uni._flipped(R) else r'position across the cell,  $z/L_z$')
        ax.set_ylabel(rf'stress  $\sigma_{{{comp}}}$  (LJ)')
        ax.grid(alpha=0.3)
        letter = f'({chr(97 + k)}) ' if n > 1 else ''
        ax.set_title(letter + (rf"free re-swelling from $\varepsilon={L['eps']:.2f}$" if ph == 'unload' else
                               rf"hold at $\varepsilon={L['eps']:.2f}$" if np.isfinite(L.get('eps', np.nan)) else 'permeation drive'), fontsize=20)
        word = 're-swelling' if ph == 'unload' else 'hold'
        handles = [Line2D([0], [0], color=SPECIES_COLORS[key], lw=3, label=rf'{nm}  $\sigma_{{{sub}{comp}}}$') for key, nm, sub in names]
        handles += [Line2D([0], [0], color='0.35', ls='--', lw=1.8, label=r'reference ($\varepsilon=0$)'),
                    Line2D([0], [0], color='0.6', lw=1.2, label=f'earlier in the {word} (pale)'),
                    Line2D([0], [0], color='0.35', lw=3, label=f'end of the {word} (bold)')]
        ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5 * sum(x_gel), 0.035), bbox_transform=tr,
                  ncol=2, fontsize=14, framealpha=0.95, borderaxespad=0.0)       # centred on the gel band, not on the axes
    fig.suptitle(f'Partial and total $\\sigma_{{{comp}}}$ profiles across the cell  |  {deck["cfg"].sim_name}', fontsize=13, fontweight='bold')
    return fig



# ---------------------------------------------------------------------------------------------
# gel volume over an equilibration run (slab_with_support deck): the postprocessing script's Rg and BB curves
# ---------------------------------------------------------------------------------------------
EQ_LOCAL = '../../flow_data_local/raw_data/chain_length'      # where the equilibration runs are synced (scripts/ relative)
VOLUME_COLORS = dict(rg='#0072B2', bb='#D55E00')              # Okabe-Ito blue / vermilion


def load_gel_volume(run, root=EQ_LOCAL, dt=0.005):
    """The two gel-volume traces of an equilibration run that scripts/postprocess.sh draws (lib/plot_lammps_log.py,
    convergence figure): gel_volume_rg (step, Rg^3 of the polymer) and gel_volume_bb (the polymer's bounding-box
    volume; no step column, same cadence), each written every 10 k steps.

    run   glob of the run directory's name under `root` (the newest match that holds a gel_volume_rg file wins);
          continuation_* folders are appended.
    Also reads the Rg-plateau bookkeeping from log.lammps when it is there (the deck's RG_HALT lines).
    -> dict(dir, step, t (tau), rg3, bb (None when the file is missing or of another length), rg_n = Rg^3/first,
            bb_n = BB/first, blocks (step, block-mean Rg) | None, detected | None, converged | None (steps)), or None."""
    import re
    import numpy as np
    root = Path(root).expanduser()

    def parts(d):
        return [d] + sorted(d.glob('continuation_*'))

    def find(d, stem):
        return [f for p in parts(d) for f in sorted(p.glob(f'output_files/volume_data/{stem}_*.dat')) + sorted(p.glob(f'{stem}_*.dat'))]

    def read(files):
        rows = []
        for f in files:
            for line in f.read_text().splitlines():
                if line and not line.startswith('#'):
                    try:
                        rows.append([float(x) for x in line.split()])
                    except ValueError:
                        pass
        return np.array(rows) if rows else np.zeros((0, 1))

    dirs = [d for d in sorted(root.glob(run)) if d.is_dir() and find(d, 'gel_volume_rg')]
    if not dirs:
        print(f"no gel_volume_rg file under {root}/{run}.  Copy the run's output_files/volume_data (and log.lammps) from the cluster, e.g.\n"
              f"  scp -r expanse:Documents/lammps_runs/<folder>/<run directory>/output_files/volume_data "
              f"'{root}/<run directory>/output_files/'")
        return None
    d = sorted(dirs, key=lambda p: p.name.rsplit('_', 2)[-2:])[-1]          # newest by the <date>_<time> stamp
    rg = read(find(d, 'gel_volume_rg'))
    bb = read(find(d, 'gel_volume_bb'))[:, 0]
    out = dict(dir=d, step=rg[:, 0], t=rg[:, 0] * dt, rg3=rg[:, 1], rg_n=rg[:, 1] / rg[0, 1], bb=None, bb_n=None,
               blocks=None, detected=None, converged=None)
    if len(bb) == len(rg):
        out.update(bb=bb, bb_n=bb / bb[0])
    elif len(bb):
        print(f'gel_volume_bb has {len(bb)} rows against {len(rg)} of gel_volume_rg: BB curve left out (it has no step column)')
    blocks = []
    for log in [p / 'log.lammps' for p in parts(d)]:
        if not log.exists():
            continue
        for line in log.read_text(errors='replace').splitlines():
            m = re.match(r'>>> RG_BLOCK step (\d+) Rg_block ([\d.]+)', line)
            if m:
                blocks.append((float(m.group(1)), float(m.group(2))))
            m = re.match(r'>>> RG_PLATEAU_DETECTED at step (\d+)', line)
            if m:
                out['detected'] = int(m.group(1))
            m = re.match(r'>>> RG_PLATEAU_CONVERGED at step (\d+)', line)
            if m:
                out['converged'] = int(m.group(1))
    out['blocks'] = np.array(blocks) if blocks else None
    return out


def fig_gel_volume(GV, layout='single', tail=0.3, plateau=True, figsize=None):
    """Gel volume over the equilibration, normalised by its first sample, as in the postprocessing script's convergence
    figure: the Rg-based volume Rg^3 and the bounding-box (BB) volume of the polymer (GV from load_gel_volume).

    layout   'single': both curves on one panel | 'stacked': one panel each on a shared time axis
    tail     fraction of the run, at its end, behind the quoted means (the postprocessing script quotes the last 30 %)
    plateau  True marks the step at which the deck's Rg-halt detected the plateau, when the log has it
    The two differ by colour (VOLUME_COLORS) and by line weight (the BB trace is the thin, noisy one)."""
    import numpy as np
    if GV is None:
        return None
    curves = [('rg', GV['rg_n'], r'$R_g^{\,3}$ volume', 2.6)]
    if GV['bb_n'] is not None:
        curves.append(('bb', GV['bb_n'], 'bounding-box (BB) volume', 1.0))
    t = GV['t'] / 1e3
    n0 = int(round(len(t) * (1.0 - tail)))
    stacked = layout == 'stacked' and len(curves) > 1
    n = len(curves) if stacked else 1
    fig, axes = plt.subplots(n, 1, figsize=figsize or (11.0, 4.3 * n + 0.8 if stacked else 6.8), sharex=True,
                             squeeze=False, constrained_layout=True)
    axes = list(axes[:, 0])
    means = []
    for k, (key, y, lab, lw) in enumerate(curves):
        ax = axes[k if stacked else 0]
        ax.plot(t, y, '-', color=VOLUME_COLORS[key], lw=lw, label=lab)
        means.append((ax, rf'{lab}: {y[n0:].mean():.4f}'))
    for ax in axes:
        if plateau and GV['detected'] is not None:
            ax.axvline(GV['detected'] * (GV['t'][0] / GV['step'][0]) / 1e3, color='0.25', ls=':', lw=1.6, label=r'$R_g$ plateau detected')
        ax.set_ylabel('gel volume / initial')
        ax.grid(alpha=0.3)
    for ax in axes:
        rows = [s for a, s in means if a is ax]
        ax.text(*((0.015, 0.06) if stacked else (0.985, 0.25)), f'mean of the last {100 * tail:.0f} %\n' + '\n'.join(rows), transform=ax.transAxes,
                ha='left' if stacked else 'right', va='bottom' if stacked else 'center', fontsize=15, bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.9))
        ax.legend(loc='upper right', fontsize=16, framealpha=0.9)
    axes[-1].set_xlabel(r'time  ($10^3\,\tau$)')
    fig.suptitle('Gel volume over the equilibration  |  ' + GV['dir'].name, fontsize=13, fontweight='bold')
    return fig
