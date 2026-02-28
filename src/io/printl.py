import datetime
from src.utils.constants import (
    VERBOSE
)

def printl(text: str):
    if VERBOSE:
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{now}] {text}")