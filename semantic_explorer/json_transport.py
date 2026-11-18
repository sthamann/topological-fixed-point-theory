def json_transport_value(value):
    """Preserve exact oversized integers as decimal strings in Atlas JSON.

    Original source data remain untouched; Rust stores these display-only
    contract attributes as generic JSON values.
    """
    if isinstance(value, dict):
        return {key: json_transport_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_transport_value(item) for item in value]
    if isinstance(value, int) and not -(1 << 63) <= value < (1 << 64):
        return str(value)
    return value
