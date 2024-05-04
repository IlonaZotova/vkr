
def lerp(x1, y1, x2, y2, x):
    return (y2 - y1) / (x2 - x1) * (x - x1) + y1
