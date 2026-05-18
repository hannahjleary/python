import matplotlib
from matplotlib import font_manager

print(matplotlib.__version__)

# for font in font_manager.findSystemFonts(fontpaths=['~/.fonts']):
#     print(font_manager.FontProperties(fname=font).get_name())

print(matplotlib.get_configdir())
print(matplotlib.get_data_path())
