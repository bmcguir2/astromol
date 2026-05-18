from astromol.database import Database

from baseline import load_production_baseline


def test_database_loads_production_data():
    db = Database()
    counts = load_production_baseline()["counts"]

    assert len(db.refs) == counts["references"]
    assert len(db.telescopes) == counts["telescopes"]
    assert len(db.sources) == counts["sources"]
    assert len(db.molecules) == counts["molecules"]
    assert len(db.detections) == counts["detections"]
