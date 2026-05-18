from astromol.database import Database


def test_database_loads_production_data():
    db = Database()

    assert len(db.refs) == 1043
    assert len(db.telescopes) == 49
    assert len(db.sources) == 92
    assert len(db.molecules) == 369
    assert len(db.detections) == 515
