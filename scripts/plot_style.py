"""Times typography for all labels and mathematical text."""
import matplotlib as mpl
from matplotlib import font_manager

def configure_times():
    for name in ('Times New Roman', 'Times', 'Nimbus Roman'):
        try:
            font_manager.findfont(font_manager.FontProperties(family=name),fallback_to_default=False)
            break
        except ValueError:
            continue
    else:
        raise RuntimeError('Install Times New Roman or the Times-compatible Nimbus Roman font.')
    mpl.rcParams.update({
        'font.family': name,
        'mathtext.fontset': 'custom',
        'mathtext.rm': name,
        'mathtext.it': name+':italic',
        'mathtext.bf': name+':bold',
        'mathtext.sf': name,
        'mathtext.tt': name,
        'mathtext.cal': name+':italic',
        'mathtext.fallback': 'stix',
    })
    return name
