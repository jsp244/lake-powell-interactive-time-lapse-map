"""Native window wrapper for the Lake Powell Timelapse viewer (pywebview/WebView2).

Loads index.html from the folder the exe lives in (falling back to the canonical
install path), so re-running fetch_imagery.py updates the app with no rebuild.
"""

import ctypes
import sys
from pathlib import Path

import webview

APP_NAME = "Lake Powell Interactive Time Lapse Map"


def base_dir():
    return Path(sys.executable).parent if getattr(sys, "frozen", False) else Path(__file__).parent


def main():
    index = base_dir() / "index.html"
    if not index.exists():
        ctypes.windll.user32.MessageBoxW(
            None,
            "index.html not found.\n\nKeep this exe in the project folder "
            "(next to index.html and frames\\).",
            APP_NAME,
            0x10,
        )
        return
    webview.create_window(
        APP_NAME,
        index.as_uri(),
        width=1280,
        height=880,
        min_size=(900, 620),
        background_color="#0d0d0d",
    )
    webview.start()


if __name__ == "__main__":
    main()
