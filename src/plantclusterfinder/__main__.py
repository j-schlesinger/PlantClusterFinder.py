"""Provide main() to provide support for `python -m plantclusterfinder`."""

from plantclusterfinder.main import main

if __name__ != "__main__":
    raise NotImplementedError("Cannot execute outside of __main__")
main()
