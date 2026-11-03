"""Validate JSON evidence at the audit's serialization boundaries."""

import json
from pathlib import Path
from typing import TypeAlias

JSONValue: TypeAlias = (
    bool | int | float | str | list["JSONValue"] | dict[str, "JSONValue"] | None
)
JSONRecord: TypeAlias = dict[str, JSONValue]


def validate_json(value: JSONValue) -> JSONValue:
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, list):
        return [validate_json(item) for item in value]
    if isinstance(value, dict):
        if not all(isinstance(key, str) for key in value):
            raise TypeError("JSON evidence requires string keys")
        return {key: validate_json(item) for key, item in value.items()}
    raise TypeError("Evidence contains a value unsupported by JSON")


def load_json(path: Path) -> JSONValue:
    return validate_json(json.loads(path.read_text()))


def record(value: JSONValue) -> JSONRecord:
    if not isinstance(value, dict):
        raise TypeError("Expected a JSON record")
    return value


def records(value: JSONValue) -> list[JSONRecord]:
    if not isinstance(value, list):
        raise TypeError("Expected a list of JSON records")
    return [record(item) for item in value]


def integer(value: JSONValue) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise TypeError("Expected an integer")
    return value


def array(value: JSONValue) -> list[JSONValue]:
    if not isinstance(value, list):
        raise TypeError("Expected a JSON array")
    return value
