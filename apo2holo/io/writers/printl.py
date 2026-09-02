import datetime
from apo2holo.utils.constants import (
    VERBOSE
)

def printl(text: str):
    if VERBOSE:
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{now}] {text}")

def print_probabilities(d, width=40):
    """
    d: dict of {label: value}
    width: max number of stars for the largest value
    """
    if VERBOSE:
        max_val = max(d.values())

        for key, val in d.items():
            bar_length = int((val / max_val) * width)
            bar = "*" * bar_length
            print(f"{str(key):<15} | {bar:<{width}} | {val:.3f}")