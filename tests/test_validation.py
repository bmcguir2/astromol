from astromol.validation import validate_database


def test_production_data_validation_has_no_errors():
    report = validate_database()

    assert report.ok
    assert report.errors == []
    assert len(report.warnings) == 5
    assert {issue.code for issue in report.warnings} == {"dipole-nonnumeric"}
    assert {issue.record for issue in report.warnings} == {
        "mol:AlCl",
        "mol:CP",
        "mol:SO+",
        "mol:MgCN",
        "mol:HNCS",
    }
