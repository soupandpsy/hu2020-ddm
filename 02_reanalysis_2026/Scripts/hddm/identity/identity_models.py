"""Upper boundary = Self. Sharing v across Self/Other avoids v/dc aliasing."""
CONDITIONS = ["val"]
UPPER = "Self"

def build(data, name, settings):
    from common import build_model
    return build_model(data, name, settings, CONDITIONS)
