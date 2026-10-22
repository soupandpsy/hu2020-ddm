"""Upper boundary = Good. Sharing v across Good/Bad avoids v/dc aliasing."""
CONDITIONS = ["id"]
UPPER = "Good"

def build(data, name, settings):
    from common import build_model
    return build_model(data, name, settings, CONDITIONS)
