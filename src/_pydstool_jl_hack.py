import scipy

original_version = scipy.__version__
try:
    scipy.__version__ = '0.9'
    import PyDSTool
finally:
    scipy.__version__ = original_version
