import runpy
import sys
import types
from pathlib import Path


def test_main_module_sets_process_title_and_exits(monkeypatch):
    calls = {}
    module_path = (
        Path(__file__).resolve().parents[1] / "src" / "anianns" / "__main__.py"
    )

    monkeypatch.setitem(
        sys.modules,
        "anianns.anianns",
        types.SimpleNamespace(main=lambda: 7),
    )
    monkeypatch.setitem(
        sys.modules,
        "setproctitle",
        types.SimpleNamespace(
            setproctitle=lambda title: calls.setdefault("title", title)
        ),
    )
    monkeypatch.setattr(
        sys, "exit", lambda code=0: (_ for _ in ()).throw(SystemExit(code))
    )

    try:
        runpy.run_path(str(module_path), run_name="__main__")
    except SystemExit as exc:
        assert exc.code == 7
    else:
        raise AssertionError("SystemExit was not raised")

    assert calls["title"] == "AniAnns"
