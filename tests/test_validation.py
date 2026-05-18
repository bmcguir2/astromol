from astromol.validation import validate_database

from baseline import load_production_baseline


def test_production_data_validation_has_no_errors():
    report = validate_database()
    baseline = load_production_baseline()
    expected_warnings = baseline["validation_warnings"]
    actual_warnings = [
        {
            "severity": issue.severity,
            "code": issue.code,
            "record": issue.record,
            "message": issue.message,
        }
        for issue in report.warnings
    ]
    actual_warnings.sort(key=lambda issue: (issue["code"], issue["record"], issue["message"]))

    assert report.ok
    assert report.errors == []
    assert actual_warnings == expected_warnings
