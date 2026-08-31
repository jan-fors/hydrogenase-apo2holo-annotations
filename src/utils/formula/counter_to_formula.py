
from collections import Counter

def counter_to_formula(counter: Counter) -> str:
    return "".join(f"{count}{key}" for key, count in sorted(counter.items()))

