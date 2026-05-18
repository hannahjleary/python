import matplotlib
from matplotlib import font_manager

font_manager.__rebuild()
print(matplotlib.get_cachedir())