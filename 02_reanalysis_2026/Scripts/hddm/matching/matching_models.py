"""Upper boundary = match choice, not Good-self. Shared v across match directions."""
CONDITIONS = ["val", "id"]
UPPER = "match"

def build(data, name, settings):
    from common import build_model
    return build_model(data, name, settings, CONDITIONS)
