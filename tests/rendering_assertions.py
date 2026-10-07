"""Presentation invariants evaluated from a test's selected data."""

def ring_text(data, *, types=False):
    categories = sorted(enumerate(data.categories), key=lambda pair: (-pair[1].count, pair[0]))
    categories = [c for _, c in categories]
    if types:
        return [c.plural_label for c in categories] + [f"{c.percent:.1f}%" for c in categories]
    return [s for c in categories for s in (c.label, f"{c.percent:.1f}%")]


def boxplot_labels(data, attr):
    keys = ("dark_cloud", "carbon_star", "sfr", "diffuse_cloud")
    categories = {c.key: c for c in data.categories}
    return [f"n={len([v for v in getattr(categories[k], attr) if attr == 'masses' or v >= 0])}"
            for k in keys if k in categories and getattr(categories[k], attr)]


def assert_slide_membership(layout, labels):
    actual = [entry.label for entry in layout.molecules]
    assert len(actual) == len(set(actual)) == layout.total
    assert set(actual) == set(labels)
    assert layout.molecule_font_pt > 0
    assert not layout.warnings
    for group in layout.groups:
        for column in group.columns:
            assert column.overflow == 0
